package srvalidate;

import beast.evolution.tree.Node;
import beast.math.statistic.DiscreteStatistics;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Alexei Drummond
 */
public class SampleTree {

    static final String S_RANGE = "stratigraphicRange";

    public SampleTree() {}

    /**
     * Takes a full tree and removes the extinct lineages
     * @param fullTree
     * @return may return a different node than the node passed if the root node itself is pruned to produce the sample tree.
     */
    public Node sampleTree(Node fullTree) {

        int m1 = markParents("reaction", "psi", fullTree);
        int m2 = markParents("rho", true, fullTree);

        //markStratigraphicIntervals(fullTree);

        markFirstAndLastOccurrences(fullTree);

//        System.out.println("nodes in full tree = " + fullTree.getAllChildNodes().size());
//        System.out.println("m1 = " + m1 + "; m2 = " + m2 + "; sum=" + (m1+m2));

        removeUnmarkedNodes(fullTree);

//        System.out.println("nodes in tree after removing unmarked nodes = " + fullTree.getAllChildNodes().size());
//
//        System.out.println("  nodes matching rho=true: " + matchesNodeTrait(fullTree.getAllChildNodes(), "rho", true));
//        System.out.println("  nodes matching rho=false: " + matchesNodeTrait(fullTree.getAllChildNodes(), "rho", false));
//        System.out.println("  nodes matching reaction=psi: " + matchesNodeTrait(fullTree.getAllChildNodes(), "reaction", "psi"));

        //fullTree = removeOlderPsiNodesAndUnobservedSpeciationNodes(fullTree);

        fullTree = removeIntermediatePsiNodesAndUnobservedSpeciationNodes(fullTree);

        removeMarkMetaData(fullTree);

        FBDTreeSimulator.processMetaData(fullTree);

        return fullTree;
    }

    private void removeMarkMetaData(Node root) {

        for (Node node : root.getAllChildNodes()) {
            node. removeMetaData("mark");
        }
    }

    /**
     *
     * @param rootNode
     * @return a new root node in the case that the original root node is removed, other returns the same root node. Either way descendant nodes removed according to rule.
     */
    private Node removeOlderPsiNodesAndUnobservedSpeciationNodes(Node rootNode) {

        List<Node> nodes = rootNode.getAllChildNodes();

        for (Node node : nodes) {
            if (isPsiSampleNode(node) && !(node.getMetaData(S_RANGE) instanceof Double) || (isSpeciation(node) && node.getChildCount()==1)) {

                Node parent = node.getParent();
                Node child = node.getChild(0);
                if (parent != null) {
                    parent.setChild(parent.getChild(0) == node ? 0 : 1, child);
                    child.setParent(parent);
                    node.setParent(null);
                } else {
                    if (node != rootNode) throw new RuntimeException("What is going on?!");
                    rootNode.removeChild(child);
                    child.setParent(null);
                    rootNode = child;
                }
            }
        }
        return rootNode;
    }

    /**
     *
     * @param rootNode
     * @return a new root node in the case that the original root node is removed, other returns the same root node. Either way descendant nodes removed according to rule.
     */
    private Node removeIntermediatePsiNodesAndUnobservedSpeciationNodes(Node rootNode) {

        List<Node> nodes = rootNode.getAllChildNodes();

        for (Node node : nodes) {
            if (isPsiSampleNode(node) && !(node.getID().contains("_first")) && !(node.getID().contains("_last")) || (isSpeciation(node) && node.getChildCount()==1)) {

                Node parent = node.getParent();
                Node child = node.getChild(0);
                if (parent != null) {
                    parent.setChild(parent.getChild(0) == node ? 0 : 1, child);
                    child.setParent(parent);
                    node.setParent(null);
                } else {
                    if (node != rootNode) throw new RuntimeException("What is going on?!");
                    rootNode.removeChild(child);
                    child.setParent(null);
                    rootNode = child;
                }
            }
        }
        return rootNode;
    }

    private void markStratigraphicIntervals(Node node) {

        // collect oldest for each rho=true
        if (node.isLeaf()) {
            Object rhoValue = node.getMetaData("rho");
            if (rhoValue instanceof Boolean &&  (Boolean)rhoValue) {
                double oldest = node.getHeight();
                Node parent = node.getParent();
                if (parent.getChild(0) == node) {
                    oldest = collectOldestSampleAge(node.getParent());
                }
                node.setMetaData(S_RANGE, Math.max(node.getHeight(),oldest)-node.getHeight());
            }
        }

        // collect oldest for each youngest psi node whose child is unmarked
        if (isYoungestPsiNodeOfSpecies(node)) {
            double oldest = node.getHeight();
            Node parent = node.getParent();
            if (parent != null && parent.getChild(0) == node) {
                oldest = collectOldestSampleAge(node.getParent());
            }
            node.setMetaData(S_RANGE, Math.max(node.getHeight(),oldest)-node.getHeight());
        }

        List<Node> children = new ArrayList<>();
        children.addAll(node.getChildren());
        for (Node child:children) {
            markStratigraphicIntervals(child);
        }
    }

    private void markFirstAndLastOccurrences(Node youngest) {

        // collect oldest for each rho=true
        if (youngest.isLeaf()) {
            Object rhoValue = youngest.getMetaData("rho");
            if (rhoValue instanceof Boolean &&  (Boolean)rhoValue) {
                String ID = youngest.getID();
                youngest.setID(ID+"_last");
                Node oldest = collectOldestSampleNode(youngest);
                if (!oldest.equals(youngest)) {
                    oldest.setID(ID+"_first");
                }

            }
        }

        // collect oldest for each youngest psi node whose child is unmarked
        if (isYoungestPsiNodeOfSpecies(youngest)) {
            String ID = youngest.getID();
            youngest.setID(ID+"_last");
            Node oldest = collectOldestSampleNode(youngest);
            if (!oldest.equals(youngest)) {
                oldest.setID(ID+"_first");
            }
        }

        List<Node> children = new ArrayList<>();
        children.addAll(youngest.getChildren());
        for (Node child:children) {
            markFirstAndLastOccurrences(child);
        }
    }

    private boolean isYoungestPsiNodeOfSpecies(Node node) {
        if (!isPsiSampleNode(node)) return false;

        Node n = node.getChild(0);
        while (isSpeciation(n)) {
            n = n.getChild(0);
        }
        return !isSampleNode(n);
    }

    public static boolean isSampleNode(Node node) {
        Object rhoValue = node.getMetaData("rho");
        return rhoValue instanceof Boolean && (Boolean) rhoValue || isPsiSampleNode(node);
    }

    public static boolean isPsiSampleNode(Node node) {
        return isReaction(node, "psi");
    }

    public static boolean isSpeciation(Node node) {
        return isReaction(node, "lambda");
    }

    public static boolean isExtinction(Node node) {
        return isReaction(node, "mu");
    }

    public static boolean isReaction(Node node, String reactionString) {
        Object reaction = node.getMetaData("reaction");
        return reaction instanceof String && ((String) reaction).equals(reactionString);
    }

    public static boolean matchesNodeTrait(Node node, String trait, Object value) {
        Object val = node.getMetaData(trait);
        return val == value;
    }

    public static int matchesNodeTrait(List<Node> nodes, String trait, Object value) {
        int matches = 0;
        for (int i = 0; i < nodes.size(); i++) {
            matches += (matchesNodeTrait(nodes.get(i), trait, value) ? 1 : 0);
        }
        return matches;
    }

    private boolean isMarked(Node node) {
        Object mark = node.getMetaData("mark");
        return mark instanceof Boolean && (Boolean) mark;
    }

    private double collectOldestSampleAge(Node node) {

        double age = 0;
        if (isSampleNode(node)) {
            age = node.getHeight();
        }

        Node parent = node.getParent();

        if (parent != null && parent.getChild(0) == node) {
            age = Math.max(age,collectOldestSampleAge(parent));
        }
        return age;
    }

    private Node collectOldestSampleNode(Node node) {
        Node oldest = null;
        if (isSampleNode(node)) {
            oldest=node;
        }

        Node parent = node.getParent();

        if (parent != null && parent.getChild(0) == node) {
            Node candidate = collectOldestSampleNode(parent);
            if (candidate != null) {
                oldest = candidate;
            }
        }

        return oldest;
    }

    private double removeOlderSamples(Node node) {

        Node parent = node.getParent();

        double age = 0;
        if (isSampleNode(node) && parent != null) {
            Node child = node.getChild(0);
            node.removeChild(child);
            parent.setChild(parent.getChild(0) == node ? 0 : 1, child);

            node.setParent(null);
            age = node.getHeight();
            // clear metadata
            node.setMetaData("reaction", null);

            node = child;
        }

        if (parent != null && parent.getChild(0) == node) {
            age = Math.max(age,removeOlderSamples(parent));
        }
        return age;
    }

    private void removeNodesByMetadata(String key, Object value, Node node) {

        List<Node> children = new ArrayList<>();
        children.addAll(node.getChildren());
        for (Node child : children) {
            removeNodesByMetadata(key, value, child);
        }

        if (node.getMetaDataNames().contains(key)) {
            Object val = node.getMetaData(key);
            if (val.equals(value)) {
                // remove this node
                Node parent = node.getParent();
                if (parent != null) {
                    parent.removeChild(node);
                }
            }
        }
    }

    private void removeUnmarkedNodes(Node node) {

        List<Node> children = new ArrayList<>();
        children.addAll(node.getChildren());
        for (Node child : children) {
            removeUnmarkedNodes(child);
        }

        if (!node.getMetaDataNames().contains("mark")) {
            // remove this node
            Node parent = node.getParent();
            if (parent != null) {
                parent.removeChild(node);
                node.setParent(null);
            }
        }
    }


    private int markParents(String key, Object value, Node node) {

        int marks = 0;

        for (Node child : node.getChildren()) {
            marks += markParents(key, value, child);
        }

        if (node.getMetaDataNames().contains(key)) {
            Object val = node.getMetaData(key);
            if (val.equals(value)) {
                // mark this node
                if (markNode(node)) marks += 1;
                Node n = node.getParent();
                while (n != null) {
                    if (markNode(n)) marks += 1;
                    n = n.getParent();
                }
            }
        }
        return marks;
    }

    /**
     * Mark the node if a mark attribute not already present
     * @param node
     * @return true if a mark attribute was added by this call, false if mark attribute already present
     */
    private boolean markNode(Node node) {
        if (node.getMetaDataNames().contains("mark")) return false;
        node.setMetaData("mark",true);
        return true;
    }

    public static void addFakeNodes(Node node) {
        if (!node.isLeaf()) {
            if(node.getChildCount() > 1) {
                for (Node child:node.getChildren()) {
                    addFakeNodes(child);
                }
            } else {
                Node child = node.getLeft();
                Node directAncestorChild = new Node();
                directAncestorChild.setID(node.getID());
                directAncestorChild.setHeight(node.getHeight());
                directAncestorChild.setParent(node);
                node.setID(null);
                node.addChild(directAncestorChild);
                addFakeNodes(child);
            }
        }
    }

    public static void main(String[] args) throws IOException {

        FileWriter fileWriterFull = new FileWriter("FullTrees.txt");
        FileWriter fileWriterSample = new FileWriter("SampleTrees.txt");
        FileWriter fileWriterTaxa = new FileWriter("SampleTaxa.txt");
        PrintWriter pwFull = new PrintWriter(fileWriterFull);
        PrintWriter pwSample = new PrintWriter(fileWriterSample);
        PrintWriter pwTaxa = new PrintWriter(fileWriterTaxa);

        double lambda = 1;
        double mu = 0.5;
        double psi = 1.5;
        double x0 = 8;
        double rho = 0.5;

        int rejectionCount = 0;
        int minTreeSize = 10;
        int reps = 1;

        double[] fullTreeSize = new double[reps];
        double[] sampleTreeSize = new double[reps];

        for (int i = 0; i < reps; i++) {

            FBDTreeSimulator simulator = new FBDTreeSimulator();

            int leaves = 0;
            Node myTree = null;

            while (leaves < minTreeSize) {
                myTree = simulator.simulate(lambda, mu, psi, x0, rho);
                leaves = myTree.getAllLeafNodes().size();
                rejectionCount += 1;
            }

            // printing full tree
            pwFull.println(myTree.toNewick() + ";");

            fullTreeSize[i] = myTree.getAllLeafNodes().size();

            SampleTree sampleTree = new SampleTree();
            myTree = sampleTree.sampleTree(myTree);

            sampleTreeSize[i] = myTree.getAllLeafNodes().size();

            addFakeNodes(myTree);

            // printing sample tree
            pwSample.println(myTree.toNewick() + ";");
            List<Node> sampledNodes = myTree.getAllLeafNodes();
            int taxonCount = sampledNodes.size();
            for (int j=0; j<taxonCount-1; j++) {
                Node sampledNode = sampledNodes.get(j);
                pwTaxa.println(sampledNode.getID()+"="+sampledNode.getHeight());
            }
            Node sampledNode = sampledNodes.get(taxonCount-1);
            pwTaxa.println(sampledNode.getID()+"="+sampledNode.getHeight() + ";");

            rejectionCount -= 1;
        }

        pwFull.flush();
        pwSample.flush();
        pwTaxa.flush();
        pwFull.close();
        pwSample.close();
        pwTaxa.close();

        System.out.println(reps + " trees simulated.");
        System.out.println(rejectionCount + " tree rejected for being too small (n<" + minTreeSize + ")");
        System.out.println("Full tree sizes (min, mean, max) = (" +
                DiscreteStatistics.min(fullTreeSize) + "," +
                DiscreteStatistics.mean(fullTreeSize) + "," +
                DiscreteStatistics.max(fullTreeSize) + ")");
        System.out.println("Sample tree sizes (min, mean, max) = (" +
                DiscreteStatistics.min(sampleTreeSize) + "," +
                DiscreteStatistics.mean(sampleTreeSize) + "," +
                DiscreteStatistics.max(sampleTreeSize) + ")");
    }

}
