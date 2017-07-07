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
     * Simulate a sample phylogeny under the stratigraphic range fossilized birth-death process with budding speciation.
     *
     * @param lambda speciation rate
     * @param mu extinction rate
     * @param psi fossil sampling rate
     * @param x0 time of the origin
     * @param rho probability of sampling an extant species at the present
     * @param unobservedSpeciationAges a list to be populated of the ages of speciation events on lineages within
     *                                  the sampled phylogeny but that are "unobserved speciation events" because only
     *                                  one of the two children lineages in the full phylogeny was subsequently sampled.
     *                                  The ages are in the same units as the node ages (i.e. time before present).
     * @param stratigraphicIntervals a list of the size of the stratigraphic intervals, indexed by node number in the
     *                               returned sample phylogeny.
     * @return a Node representing the root of the sample phylogeny after extinct lineages have been pruned away.
     */
    public Node simulate(double lambda, double mu, double psi, double x0, double rho,
                         List<Double> unobservedSpeciationAges, List<Double> stratigraphicIntervals) {


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
                youngestFossil.setMetaData("oldestFossil",(Double)srange.getValue(0));
                youngestFossil.setMetaData("youngestFossil",(Double)srange.getValue(1));
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
                node.setMetaData("oldestFossil", (Double)srange.getValue(0));
                node.setMetaData("youngestFossil", (Double)srange.getValue(1));
            }
        }

        // remove metadata from unobserved speciation and remove the nodes from the srange map
        for (Node node : toBeRemovedFromSranges) {
            node.removeMetaData("oldestFossil");
            node.removeMetaData("youngestFossil");
            srangeMap.remove(node);
        }

        // remove the unobserved speciation nodes in the sample tree and store their ages in the provided list.
        removeUnobservedSpeciationNodes(start, unobservedSpeciationAges);

        processMetaData(start);

        int numSRanges = srangeMap.values().size();

        int i = 0;
        int j = numSRanges;
        for (Node node : start.getAllChildNodes()) {
            if (isStratigraphicRangeEnd(node,srangeMap)) {
                node.setNr(i);
                node.setID(i+"");
                RealParameter srange = srangeMap.get(node);
                stratigraphicIntervals.add(srange.getValue(0) - srange.getValue(1));
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
     * @param ages a list of node heights to accumulate, one for each unobserved speciation node removed.
     */
    private void removeUnobservedSpeciationNodes(Node node, List<Double> ages) {
        List<Node> children = new ArrayList<>();
        children.addAll(node.getChildren());
        for (Node child : children) {
            removeUnobservedSpeciationNodes(child,ages);
        }
        if (isUnobservedSpeciationNode(node)) {
            // remove it and store times
            ages.add(node.getHeight());

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
        return node.getParent() != null && node.getChildCount() == 1 &&
                (getDoubleMetaData(node, "youngestFossil") == null || getDoubleMetaData(node, "youngestFossil") < node.getHeight());
    }

    private Double getDoubleMetaData(Node node, String key) {

        Object o = node.getMetaData(key);
        if (o instanceof Integer && (Integer)o == 0) {
            return null;
        }
        return (Double)o;
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

        int n = 0;

        List<Double> unobservedSpeciationAges = new ArrayList<>();
        List<Double> stratigraphicIntervals = new ArrayList<>();
        Node myTree = null;
        while (n < 10) {
            unobservedSpeciationAges.clear();
            stratigraphicIntervals.clear();
            myTree = simulator.simulate(lambda, mu, psi, x0, rho, unobservedSpeciationAges, stratigraphicIntervals);
            n = myTree.getAllLeafNodes().size();
        }

        // clear meta data strings and remove speciation node ids
        clearMetaDataStringsRemoveSpeciationNodeIds(myTree);

        writeDensityMapperXML("myxml.xml", lambda, mu, psi, x0, rho, myTree, stratigraphicIntervals, unobservedSpeciationAges);
    }

    private static void clearMetaDataStringsRemoveSpeciationNodeIds(Node node) {
        for (Node child : node.getChildren()) {
            clearMetaDataStringsRemoveSpeciationNodeIds(child);
        }
        node.metaDataString = null;
        // if this is a speciation node then remove the id.
        if (node.getChildCount() == 2) {
            node.setID(null);
        }
    }

    private static void writeDensityMapperXML(String xmlFileName,
                                              double lambda,
                                              double mu,
                                              double psi,
                                              double x0,
                                              double rho,
                                              Node myTree,
                                              List<Double> stratigraphicIntervals,
                                              List<Double> unobservedSpeciationAges) throws IOException {

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
                "        <realParam spec=\"RealParameter\" id=\"psi\" value=\"" + psi +
                "\" lower=\"0.1\" upper=\"3.0\"/>\n" +
                "        <steps spec=\"IntegerParameter\" value=\"80\"/>\n" +
                "\n" +
                "        <distribution spec=\"SRTreeDensity\" id=\"density\">\n");

        writer.write("            <lambda spec=\"RealParameter\" value=\"" + lambda + "\"/>\n");
        writer.write("            <mu spec=\"RealParameter\" value=\"" + mu + "\"/>\n");
        writer.write("            <psi idref=\"psi\"/>");
        writer.write("            <x0 spec=\"RealParameter\" value=\"" + x0 + "\"/>\n");
        writer.write("            <rho spec=\"RealParameter\" value=\"" + rho + "\"/>\n");

        writer.write("            <tree spec=\"TreeParser\" adjustTipHeights=\"false\" newick=\"" + myTree.toNewick()+
                "\" offset=\"0\"/>\n");
        writer.write("            <sranges spec=\"RealParameter\" dimension=\"" + stratigraphicIntervals.size() +
                "\" value=\"" + spaceDelimited(stratigraphicIntervals) + "\"/>\n");

        writer.write("            <unobsSpecTimes spec=\"RealParameter\" dimension=\"" + unobservedSpeciationAges.size()
                + "\" value=\"" + spaceDelimited(unobservedSpeciationAges) + "\"/>\n");

        writer.write("        </distribution>\n" +
                "\n" +
                "        <logger fileName=\"$(filebase).log\" logEvery=\"1\">\n" +
                "            <log idref=\"psi\"/>\n" +
                "            <log idref=\"density\"/>\n" +
                "        </logger>\n" +
                "\n" +
                "        <logger id=\"screenlog\" logEvery=\"1\">\n" +
                "            <log idref=\"psi\"/>\n" +
                "            <log idref=\"density\"/>\n" +
                "        </logger>\n" +
                "    </run>\n" +
                "</beast>\n");

        writer.flush();
        writer.close();
    }

    private static String spaceDelimited(List<Double> doubleList) {
        StringBuilder builder = new StringBuilder();
        builder.append(doubleList.get(0));
        for (int i = 1; i < doubleList.size(); i++) {
            builder.append(" ");
            builder.append(doubleList.get(i));
        }
        return builder.toString();
    }
}
