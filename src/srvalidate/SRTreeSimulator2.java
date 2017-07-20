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
 * <p>
 * Mathematical model details can be found:
 * <p>
 * "The fossilized birth-death model for the analysis of stratigraphic range data under different speciation concepts" (2017)
 * by Tanja Stadler, Alexandra Gavryushkina, Rachel C.M. Warnock, Alexei J. Drummond, Tracy A. Heath
 * https://arxiv.org/abs/1706.10106
 * <p>
 * The model implemented here is the budding speciation model of theorem 7, producing oriented sampled trees.
 *
 * @author Alexei Drummond
 */
public class SRTreeSimulator2 {


    /**
     * Simulate a sample phylogeny under the stratigraphic range fossilized birth-death process with budding speciation.
     *
     * @param lambda                   speciation rate
     * @param mu                       extinction rate
     * @param psi                      fossil sampling rate
     * @param x0                       time of the origin
     * @param rho                      probability of sampling an extant species at the present
     * @param unobservedSpeciationAges a list to be populated of the ages of speciation events on lineages within
     *                                 the sampled phylogeny but that are "unobserved speciation events" because only
     *                                 one of the two children lineages in the full phylogeny was subsequently sampled.
     *                                 The ages are in the same units as the node ages (i.e. time before present).
     * @param stratigraphicIntervals   a list of the size of the stratigraphic intervals, indexed by node number in the
     *                                 returned sample phylogeny.
     * @return a Node representing the root of the sample phylogeny after extinct lineages have been pruned away.
     */
    public Node simulate(double lambda, double mu, double psi, double x0, double rho,
                         List<Double> unobservedSpeciationAges, List<Double> stratigraphicIntervals) {


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

            if (U < totalBirthRate / totalPropensity) {
                // do birth
                doBirth(activeLineages, T, x0, taxa);

            } else {
                U -= totalBirthRate / totalPropensity;
                if (U < totalDeathRate / totalPropensity) {
                    // do death
                    doDeath(activeLineages, T, x0);

                } else {
                    // do sample
                    Node sampledNode = activeLineages.get(Randomizer.nextInt(activeLineages.size()));
                    doSample(sampledNode, T, x0, taxa);
                }
            }
        }

        for (Node activeLineage : activeLineages) {
            if (Randomizer.nextDouble() < rho) {
                doSample(activeLineage, x0, x0, taxa);
            }
        }

        // At this point in the process we have produced a full tree, containing both sampled and unsampled lineages.

        // traverse the whole tree, removing subtrees that have no stratigraphic intervals in them
        int nodesRemoved = removeUnsampled(start);

        List<Node> youngestSampledNodes = new ArrayList<>();
        getAllYoungestSampledNodes(start, youngestSampledNodes);

        int i = 0;
        for (Node node : youngestSampledNodes) {
            node.setNr(i);
            node.setID(i + "");
            double y = node.getHeight();
            double o = getSRangeEnd(node);
            stratigraphicIntervals.add(o - y);
            i += 1;
        }

        // remove the unobserved speciation nodes in the sample tree and store their ages in the provided list.
        removeUnobservedSpeciationNodes(start, unobservedSpeciationAges);

        cleanup(start);

        return start;
    }

    private void getAllYoungestSampledNodes(Node node, List<Node> youngestNodes) {

        if (isSampled(node)) {
            if (node.getChildCount() == 0) youngestNodes.add(node);

            if (node.getChildCount() == 1 && isUnobservedSpeciationNode(node.getChild(0))) youngestNodes.add(node);
        }

        for (Node child : node.getChildren()) {
            getAllYoungestSampledNodes(child, youngestNodes);
        }
    }

    /**
     * @param node
     * @return the age of the oldest sampled node on the left branch above this sampled node
     */
    private double getSRangeEnd(Node node) {
        double o = node.getHeight();
        Node n = node;
        while (isLeft(n) && n.getParent() != null) {
            n = n.getParent();
            if (isSampled(n)) {
                o = n.getHeight();
            }
        }
        return o;
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
            removeUnobservedSpeciationNodes(child, ages);
        }
        if (isUnobservedSpeciationNode(node)) {
            // remove it and store times
            ages.add(node.getHeight());

            Node parent = node.getParent();

            if (parent != null) {
                Node child = node.getChild(0);
                parent.removeChild(node);
                node.removeChild(child);
                parent.addChild(child);
            }
        }
    }

    /**
     * @param node the node to test
     * @return true if the node is an unobserved speciation node.
     */
    private boolean isUnobservedSpeciationNode(Node node) {
        boolean isUnobserved1 = (node.getParent() != null && node.getChildCount() == 1 && !isSampled(node));
        boolean isUnobserved2 = node.getMetaData("childRemoved") instanceof String &&
                (node.getMetaData("childRemoved").equals("left") || node.getMetaData("childRemoved").equals("right"));
//        if (isUnobserved1 != isUnobserved2) {
//            System.out.println("node.getParent()="+node.getParent());
//            System.out.println("node.getChildCount()=" + node.getChildCount());
//            System.out.println("isSampled(node)="+isSampled(node));
//            System.out.println("node.childRemoved="+node.getMetaData("childRemoved"));
//            throw new RuntimeException("This shouldn't happen!");
//        }
        return isUnobserved2;
    }

    /**
     * This method recursively removes all "unsampled nodes" in a post-order traversal.
     *
     * @param node the node to consider for removal, after having recursively considered it children for removal.
     * @return the number of unsampled nodes removed
     */
    private int removeUnsampled(Node node) {

        int removed = 0;
        List<Node> children = new ArrayList<>();
        children.addAll(node.getChildren());

        for (Node child : children) {
            removed += removeUnsampled(child);
        }

        if (node.isLeaf() && !isSampled(node)) {

            if (node.getParent() != null) {
                Node parent = node.getParent();
                parent.setMetaData("childRemoved", isLeft(node) ? "left" : "right");
                parent.removeChild(node);
                removed += 1;
            }
        }
        return removed;
    }

    private boolean isSampled(Node node) {
        Object isSampled = node.getMetaData("sampled");
        if (isSampled instanceof Integer && (Integer)isSampled == 0) return false;
        return (boolean)isSampled;
    }

    /**
     * Creates a new node. Adds the unique id for this node to the given set of taxa.
     *
     * @param taxa a set of taxa that already exists in previously created nodes.
     * @return a new node with a new unique taxon with node number = taxa.size()+1.
     */
    private Node createNode(Set<String> taxa) {
        String newTaxon = (taxa.size() + 1) + "";
        if (taxa.contains(newTaxon)) throw new RuntimeException("Expecting taxon " + newTaxon + " to be new!");
        taxa.add(newTaxon);
        Node newNode = new Node(newTaxon);
        newNode.setNr(taxa.size() - 1);
        return newNode;
    }

    /**
     * Perform a speciation event. This creates a new node at the current time and adds it to the active lineages list.
     *
     * @param nodes the set of nodes representing "active lineages"
     * @param time  the current simulation time (since the time of origin).
     * @param x0    the age of the origin before the present.
     * @param taxa  a set of taxa that already exist.
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
     *
     * @param nodes the set of nodes representing active lineages.
     */
    private void doDeath(List<Node> nodes, double time, double x0) {

        Node deadNode = nodes.get(Randomizer.nextInt(nodes.size()));
        nodes.remove(deadNode);
    }

    /**
     * Returns true if the node is the logically left child or the root.
     * For children of unobserved speciation nodes, checks if parent.childRemoved == "right";
     * Returns false if the node is null.
     *
     * @param node the node to test.
     * @return true if the node is the left child or the root. Returns false if the node is null.
     */
    private static boolean isLeft(Node node) {
        if (node == null) return false;
        if (node.getParent() == null) return true;

        Node parent = node.getParent();
        Object childRemoved = parent.getMetaData("childRemoved");
        if (childRemoved != null) {
            return childRemoved.equals("right");
        }

        return parent.getChild(0) == node;
    }

    /**
     * Perform a fossil sampling event. This samples the given lineage by creating a new child node and
     * associating it with a new or existing stratigraphic range depending on whether the parent already represents
     * a stratigraphic range branch.
     *
     * @param nodeToSample the lineage to create a fossil sample from.
     * @param time         the current simulation time (since the time of origin).
     * @param x0           the age of the origin before the present.
     */
    private void doSample(Node nodeToSample, double time, double x0, Set<String> taxa) {

        Node sampleNode = createNode(taxa);

        sampleNode.setHeight(x0 - time);
        sampleNode.setMetaData("sample", true);

        nodeToSample.addChild(sampleNode);
    }

    public static void main(String[] args) throws IOException {

        SRTreeSimulator2 simulator = new SRTreeSimulator2();

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

        // clear meta data strings, remove speciation node ids and add fake nodes

        writeDensityMapperXML("myxml.xml", lambda, mu, psi, x0, rho, myTree, stratigraphicIntervals, unobservedSpeciationAges);

        // then run myxml.xml in BeastMain and plot in R:
        // library(ggplot2)
        // val <- read.table("~/Git/sranges_validation/examples/myxml.log", sep="\t", header=T)
        // ggplot(val, aes(lambda, psi)) + geom_raster(aes(fill = exp(density+max(density))))
    }

    private static void cleanup(Node node) {
        for (Node child : node.getChildren()) {
            cleanup(child);
        }
        node.metaDataString = null;
        // if this is sampled ancestor add a fake node
        if (node.getChildCount() == 1) {
            Node parent = node.getParent();
            Node child = node.getChild(0);

            Node fake = new Node();
            fake.setHeight(node.getHeight());

            node.removeChild(child);

            if (parent != null) {
                parent.setChild(isLeft(node) ? 0 : 1, fake);
            }
            fake.addChild(child);
            fake.addChild(node);
        } else if (node.getChildCount() == 2) {
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
                "        <realParam spec=\"RealParameter\" id=\"lambda\" value=\"" + lambda +
                "\" lower=\"0.1\" upper=\"3.0\"/>\n" +
                "        <realParam spec=\"RealParameter\" id=\"psi\" value=\"" + psi +
                "\" lower=\"0.1\" upper=\"3.0\"/>\n" +
                "        <steps spec=\"IntegerParameter\" value=\"80\"/>\n" +
                "        <steps spec=\"IntegerParameter\" value=\"80\"/>\n" +
                "\n" +
                "        <distribution spec=\"SRTreeDensity\" id=\"density\">\n");

        writer.write("            <lambda idref=\"lambda\"/>\n");
        writer.write("            <mu spec=\"RealParameter\" value=\"" + mu + "\"/>\n");
        writer.write("            <psi idref=\"psi\"/>\n");
        writer.write("            <x0 spec=\"RealParameter\" value=\"" + x0 + "\"/>\n");
        writer.write("            <rho spec=\"RealParameter\" value=\"" + rho + "\"/>\n");

        writer.write("            <tree spec=\"TreeParser\" adjustTipHeights=\"false\" IsLabelledNewick=\"true\" newick=\"" + myTree.toNewick() +
                "\" offset=\"0\"/>\n");
        writer.write("            <sranges spec=\"RealParameter\" dimension=\"" + stratigraphicIntervals.size() +
                "\" value=\"" + spaceDelimited(stratigraphicIntervals) + "\"/>\n");

        writer.write("            <unobsSpecTimes spec=\"RealParameter\" dimension=\"" + unobservedSpeciationAges.size()
                + "\" value=\"" + spaceDelimited(unobservedSpeciationAges) + "\"/>\n");

        writer.write("        </distribution>\n" +
                "\n" +
                "        <logger fileName=\"$(filebase).log\" logEvery=\"1\">\n" +
                "            <log idref=\"psi\"/>\n" +
                "            <log idref=\"lambda\"/>\n" +
                "            <log idref=\"density\"/>\n" +
                "        </logger>\n" +
                "\n" +
                "        <logger id=\"screenlog\" logEvery=\"1\">\n" +
                "            <log idref=\"psi\"/>\n" +
                "            <log idref=\"lambda\"/>\n" +
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