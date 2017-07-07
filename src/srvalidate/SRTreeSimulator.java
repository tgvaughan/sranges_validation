package srvalidate;

import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.util.Randomizer;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

/**
 * A class that simulates trees from the stratigraphic range fossilized birth-death process.
 *
 * Mathematical model details can be found:
 *
 * "The fossilized birth-death model for the analysis of stratigraphic range data under different speciation concepts" (2017)
 * by Tanja Stadler, Alexandra Gavryushkina, Rachel C.M. Warnock, Alexei J. Drummond, Tracy A. Heath
 * https://arxiv.org/abs/1706.10106
 *
 * The model implemented here is the budding speciation model of theorem 7, producing oriented sampled trees.
 *
 * @author Alexei Drummond
 * @author Laura Garcia
 */
public class SRTreeSimulator {

    /**
     *
     * @param lambda speciation rate
     * @param mu extinction rate
     * @param psi fossil sampling rate
     * @param x0 time of the origin
     * @param rho probability of sampling an extant species at the present
     * @return a Node representing the root of the sample phylogeny after extinct lineages have been pruned away.
     */
    public Node simulate(double lambda, double mu, double psi, double x0, double rho) {


        Map<Node, RealParameter> srangeMap = new HashMap<>();
        Set<String> taxa = new TreeSet<>();

        double T = 0;

        Node start = createNode(taxa);
        start.setHeight(x0);

        List<Node> activeLineages = new ArrayList<>();

        activeLineages.add(start);

        while (T < x0) {

            int n = activeLineages.size();

            double totalBirthRate = lambda * n;
            double totalDeathRate = mu * n;
            double totalSampleRate = psi * n;

            double totalPropensity = totalBirthRate + totalDeathRate + totalSampleRate;

            double dt = Randomizer.nextExponential(totalPropensity);

            double newT = T + dt;

            if (newT > x0) {
                T = x0;
                break;
            }

            T = newT;
            double U = Randomizer.nextDouble();

            if (U < totalBirthRate/totalPropensity) {
                // do birth
                doBirth(activeLineages, T, x0, taxa);

            } else {
                U -= totalBirthRate/totalPropensity;
                if (U < totalDeathRate/totalPropensity) {
                    // do death
                    doDeath(activeLineages, T, x0);

                } else {
                    // do sample
                    Node sampledNode = activeLineages.get(Randomizer.nextInt(activeLineages.size()));
                    doSample(sampledNode, T, x0, srangeMap);
                }
            }
        }

        for (Node activeLineage : activeLineages) {
            if (Randomizer.nextDouble() < rho) {
                doSample(activeLineage, x0, x0, srangeMap);
            }
        }

        // At this point in the process we have produced a full tree, containing both sampled and unsampled lineages.

        double epsilon = 1e-8;

        // mark all nodes that are part of stratigraphic ranges with their oldest and youngest fossil ages metadata
        for (Node node : srangeMap.keySet()) {
            RealParameter srange = srangeMap.get(node);
            node.setMetaData("oldestFossil",srange.getValue(0));
            node.setMetaData("youngestFossil",srange.getValue(1));
        }

        // traverse the whole tree, removing subtrees that have no stratigraphic intervals in them
        int nodesRemoved = removeUnsampled(start, srangeMap);


        //set leaves to the height of the youngest fossil in their stratigraphic range
        // this is necessary because an internal node with two extinct lineages for children may become a leaf
        // after trimming extinct lineages that needs its age needs to be increased to the youngest fossil of the
        // stratigraphic range associated with it.
        for (Node leaf : start.getAllLeafNodes()) {
            leaf.setHeight(srangeMap.get(leaf).getValue(1));
        }


        List<Node> toBeRemovedFromSranges = new ArrayList<>();
        Set<Node> keys = new HashSet<>();
        keys.addAll(srangeMap.keySet());

        // above each unobserved speciation node closest to a stratigraphic range above,
        // introduce a new node representing the youngest fossil in the stratigraphic range.
        for (Node node : keys) {
            RealParameter srange = srangeMap.get(node);

            if (!isStratigraphicRangeEnd(node,srangeMap) && node.getHeight() < srange.getValue(1)) {
                // introduce new node for youngest fossil
                Node youngestFossil = createNode(taxa);
                youngestFossil.setHeight(srange.getValue(1));
                youngestFossil.setMetaData("oldestFossil",srange.getValue(0));
                youngestFossil.setMetaData("youngestFossil",srange.getValue(1));
                Node parent = node.getParent();
                if (parent != null) {
                    parent.removeChild(node);
                    parent.addChild(youngestFossil);
                    youngestFossil.addChild(node);
                } else if (node == start) {
                    youngestFossil.addChild(node);
                    start = youngestFossil;
                } else throw new RuntimeException("No idea what happened here!");

                toBeRemovedFromSranges.add(node);
                // add in new node to map
                srangeMap.put(youngestFossil, srange);

            } else {
                node.setMetaData("oldestFossil", srange.getValue(0));
                node.setMetaData("youngestFossil", srange.getValue(1));
            }
        }

        // remove metadata from unobserved speciation and remove the nodes from the srange map
        for (Node node : toBeRemovedFromSranges) {
            node.removeMetaData("oldestFossil");
            node.removeMetaData("youngestFossil");
            srangeMap.remove(node);
        }

        List<Double> speciationTimes = new ArrayList<>();
        removeUnobservedSpeciationNodes(start, speciationTimes);
        System.out.println("unobserved speciation times:"+speciationTimes);

        processMetaData(start);

        int numSRanges = srangeMap.values().size();

        int i = 0;
        int j = numSRanges;
        for (Node node : start.getAllChildNodes()) {
            if (isStratigraphicRangeEnd(node,srangeMap)) {
                node.setNr(i);
                node.setID(i+"");
                RealParameter srange = srangeMap.get(node);
                System.out.print(" "+(srange.getValue(0)-srange.getValue(1)));
                i += 1;
            } else {
                node.setNr(j);
                node.setID(j+"");
                j += 1;
            }
        }

        return start;
    }

    /**
     * This method recursively removes all "unobserved speciation nodes" in a post-order traversal.
     *
     * @param node the node to consider for removal, after having recursively considered its children for removal.
     * @param times a list of node times to accumulate, one for each node removed.
     */
    private void removeUnobservedSpeciationNodes(Node node, List<Double> times) {
        List<Node> children = new ArrayList<>();
        children.addAll(node.getChildren());
        for (Node child : children) {
            removeUnobservedSpeciationNodes(child,times);
        }
        if (isUnobservedSpeciationNode(node)) {
            // remove it and store times
            times.add(node.getHeight());

            Node parent = node.getParent();
            Node child = node.getChild(0);

            parent.removeChild(node);
            node.removeChild(child);
            parent.addChild(child);
        }
    }

    /**
     * @param node the node to test
     * @return true if the node is an unobserved speciation node.
     */
    private boolean isUnobservedSpeciationNode(Node node) {
        return node.getParent() != null && node.getChildCount() == 1 && (node.getMetaData("youngestFossil") == null || (Double)node.getMetaData("youngestFossil") < node.getHeight());
    }

    /**
     *
     * @param node the node to test
     * @param srangeMap the map from nodes to stratigraphic ranges
     * @return true if this node represents the most recent fossil observation of a stratigraphic range (i.e. it is a y_i node in the paper).
     */
    private boolean isStratigraphicRangeEnd(Node node, Map<Node, RealParameter> srangeMap) {
        RealParameter srange = srangeMap.get(node);
        if (srange == null) return false;
        return Math.abs(srange.getValue(1) - node.getHeight()) < 1e-8;
    }

    /**
     * This method recursively removes all "unsampled nodes" in a post-order traversal.
     *
     * @param node the node to consider for removal, after having recursively considered it children for removal.
     * @param srangeMap the map from nodes to stratigraphic ranges
     * @return the number of unsampled nodes removed
     */
    private int removeUnsampled(Node node, Map<Node, RealParameter> srangeMap) {

        int removed = 0;
        List<Node> children = new ArrayList<>();
        children.addAll(node.getChildren());

        for (Node child : children) {
            removed += removeUnsampled(child,srangeMap);
        }

        if (node.isLeaf() && srangeMap.get(node) == null) {

            if (node.getParent() != null) {
                Node parent = node.getParent();
                parent.removeChild(node);
                parent.setMetaData("childRemoved", true);
                removed += 1;
            }
        }
        return removed;
    }

    /**
     * Creates a new node. Adds the unique id for this node to the given set of taxa.
     *
     * @param taxa a set of taxa that already exists in previously created nodes.
     * @return a new node with a new unique taxon with node number = taxa.size()+1.
     */
    private Node createNode(Set<String> taxa) {
        String newTaxon = (taxa.size()+1)+"";
        if (taxa.contains(newTaxon)) throw new RuntimeException("Expecting taxon "+newTaxon+" to be new!");
        taxa.add(newTaxon);
        Node newNode = new Node(newTaxon);
        newNode.setNr(taxa.size()-1);
        return newNode;
    }

    /**
     * Perform a speciation event. This creates a new node at the current time and adds it to the active lineages list.
     *
     * @param nodes the set of nodes representing "active lineages"
     * @param time the current simulation time (since the time of origin).
     * @param x0 the age of the origin before the present.
     * @param taxa a set of taxa that already exist.
     */
    private void doBirth(List<Node> nodes, double time, double x0, Set<String> taxa) {

        Node parent = nodes.get(Randomizer.nextInt(nodes.size()));
        parent.setHeight(x0 - time);

        Node leftChild = createNode(taxa);
        Node rightChild = createNode(taxa);

        leftChild.setHeight(x0 - time);
        rightChild.setHeight(x0 - time);

        parent.addChild(leftChild);
        parent.addChild(rightChild);

        nodes.remove(parent);

        nodes.add(leftChild);
        nodes.add(rightChild);
    }

    /**
     * Performs an extinction event. This pick a random active lineage and removes it from the active lineage list.
     * @param nodes the set of nodes representing active lineages.
     */
    private void doDeath(List<Node> nodes, double time, double x0) {

        Node deadNode = nodes.get(Randomizer.nextInt(nodes.size()));
        nodes.remove(deadNode);

        // should I set the nodes time to time
        // deadNode.setHeight(x0 - time);
    }

    /**
     * Returns true if the node is the left child or the root. Returns false if the node is null.
     * @param node the node to test.
     * @return true if the node is the left child or the root. Returns false if the node is null.
     */
    private boolean isLeft(Node node) {
        if (node == null) return false;
        if (node.getParent() == null) return true;
        return node.getParent().getChild(0) == node;
    }

    /**
     * Perform a fossil sampling event. This samples the given lineage by creating a new child node and
     * associating with a new or existing stratigraphic range depending on whether the parent already represents
     * a stratigraphic range branch.
     * @param nodeToSample the lineage to create a fossil sample from.
     * @param time the current simulation time (since the time of origin).
     * @param x0 the age of the origin before the present.
     * @param srangeMap the map of nodes to stratigraphic intervals.
     */
    private void doSample(Node nodeToSample, double time, double x0, Map<Node,RealParameter> srangeMap) {

        RealParameter srange = srangeMap.get(nodeToSample);
        if (srange == null && isLeft(nodeToSample)) {
            Node parent = nodeToSample.getParent();
            while (parent != null) {
                srange = srangeMap.get(parent);
                if (srange != null) break;
                if (isLeft(parent)) {
                    parent = parent.getParent();
                } else parent = null;
            }
            if (srange != null) {
                Node node = parent;
                do {
                    node = node.getChild(0);
                    srangeMap.put(node, srange);
                } while (node != nodeToSample);
            }
        }

        if (srange == null) {
            RealParameter newsrange = new RealParameter(new Double[2]);
            newsrange.setValue(0,x0-time);
            newsrange.setValue(1,x0-time);
            srangeMap.put(nodeToSample,newsrange);
        } else {
            srange.setValue(1,x0-time);
        }
        nodeToSample.setHeight(x0-time);
    }

    /**
     * This method recursively converts each node's metadata object to the metadata string. Ideally this should
     * be implemented in the core.
     * @param node the node to process the metadata of.
     */
    private void processMetaData(Node node) {
        for (Node child : node.getChildren()) {
            processMetaData(child);
        }
        Set<String> metaDataNames = node.getMetaDataNames();
        if (metaDataNames != null && !metaDataNames.isEmpty()) {
            String metadata = "";
            for (String name : metaDataNames) {
                Object value = node.getMetaData(name);
                metadata += name + "=";
                if (value instanceof Object[]) {
                    Object [] values = (Object[]) value;
                    metadata += "{";
                    for (int i = 0; i < values.length; i++) {
                        metadata += values[i].toString();
                        if (i < values.length - 1) {
                            metadata += ",";
                        }
                    }
                    metadata += "}";
                } else {
                    metadata += value.toString();
                }
                metadata += ",";
            }
            metadata = metadata.substring(0, metadata.length() - 1);
            node.metaDataString = metadata;
        }
    }


    public static void main(String[] args) throws IOException {

        SRTreeSimulator simulator = new SRTreeSimulator();

        double lambda = 1;
        double mu = 0.5;
        double psi = 2.0;
        double x0 = 8;
        double rho = 0.5;

//        String fileName = "SRTreeSimulator_l"+lambda+"_mu" + mu + "_psi" + psi + "_T" + x0 + "_p" + rho + ".txt";

        //       File outFile = new File(fileName);
//        PrintWriter writer = new PrintWriter(new FileWriter(outFile));

//        writer.println("rep\tnumLeaves");
        ////       int count = 0;
        ////       int reps = 10;
        ////       int minSize = 10;

        ////      while (count < reps) {
        Node myTree = simulator.simulate(lambda, mu, psi, x0, rho);

        String treeInNewick = myTree.toNewick();
        System.out.println(treeInNewick);

        ////          if (myTree.getAllLeafNodes().size() > minSize) {
        ////              String xmlFileName = "SRTreeDensity_l"+lambda+"_mu" + mu + "_psi" + psi + "_T" + x0 + "_p" + rho + "_" + count + ".xml";

        ////              writeSRTreeDensityXML(xmlFileName, lambda, mu, psi, x0, rho, myTree);
        ////              count += 1;
        ////          }

        ////}
//        writer.flush();
//        writer.close();
    }

    private static void writeSRTreeDensityXML(String xmlFileName, double lambda, double mu, double psi, double x0, double rho, Node myTree) throws IOException {

        File outFile = new File(xmlFileName);
        PrintWriter writer = new PrintWriter(new FileWriter(outFile));

        writer.write("<beast version=\"2.0\"\n" +
                "    namespace=\"\n" +
                "        beast.util\n" +
                "        :beast.core.parameter\n" +
                "        :beast.evolution.tree.coalescent\n" +
                "        :feast.mapping\n" +
                "        :srvalidate\">\n" +
                "    <run spec=\"DensityMapper\">\n" +
                "\n" +
                "        <realParam spec=\"RealParameter\" id=\"x0\" value=\"0.5\" lower=\"0.1\" upper=\"8.0\"/>\n" +
                "        <steps spec=\"IntegerParameter\" value=\"80\"/>\n" +
                "\n" +
                "        <distribution spec=\"SRTreeDensity\" id=\"density\">\n");

        writer.write("            <lambda spec=\"RealParameter\" value=\"" + lambda + "\"/>\n");

        // TODO the rest :)

        writer.flush();
        writer.close();
    }
}
