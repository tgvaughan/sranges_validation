package srvalidate;

import beast.evolution.tree.Node;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Alexei Drummond
 */
public class SampleTree {

    public SampleTree() {}

    /**
     * Takes a full tree and removes the extinct lineages
     * @param fullTree
     */
    public void sampleTree(Node fullTree) {

        markParents("reaction", "psi", fullTree);
        markParents("rho", true, fullTree);

        markStratigraphicIntervals(fullTree);

        removeUnmarkedNodes(fullTree);

        removeOlderPsiNodesAndUnobservedSpeciationNodes(fullTree);

        FBDTreeSimulator.processMetaData(fullTree);
    }

    private void removeOlderPsiNodesAndUnobservedSpeciationNodes(Node rootNode) {

        List<Node> nodes = rootNode.getAllChildNodes();

        for (Node node : nodes) {
            if (isPsiSampleNode(node) && !(node.getMetaData("oldestFossil") instanceof Double) || (isSpeciation(node) && node.getChildCount()==1)) {

                Node parent = node.getParent();
                Node child = node.getChild(0);
                if (parent != null) {
                    parent.setChild(parent.getChild(0) == node ? 0 : 1, child);
                    child.setParent(parent);
                    node.setParent(null);
                }
            }
        }
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
                node.setMetaData("oldestFossil", Math.max(node.getHeight(),oldest));
            }
        }

        // collect oldest for each youngest psi node whose child is unmarked
        if (isYoungestPsiNodeOfSpecies(node)) {
            double oldest = node.getHeight();
            Node parent = node.getParent();
            if (parent.getChild(0) == node) {
                oldest = collectOldestSampleAge(node.getParent());
            }
            node.setMetaData("oldestFossil", Math.max(node.getHeight(),oldest));
        }

        List<Node> children = new ArrayList<>();
        children.addAll(node.getChildren());
        for (Node child:children) {
            markStratigraphicIntervals(child);
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

    private boolean isSampleNode(Node node) {
        Object rhoValue = node.getMetaData("rho");
        return rhoValue instanceof Boolean && (Boolean) rhoValue || isPsiSampleNode(node);
    }

    private boolean isPsiSampleNode(Node node) {
        return isReaction(node, "psi");
    }

    private boolean isSpeciation(Node node) {
        return isReaction(node, "lambda");
    }

    private boolean isExtinction(Node node) {
        return isReaction(node, "mu");
    }

    private boolean isReaction(Node node, String reactionString) {
        Object reaction = node.getMetaData("reaction");
        return reaction instanceof String && ((String) reaction).equals(reactionString);
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
            }
        }
    }


    private void markParents(String key, Object value, Node node) {

        for (Node child : node.getChildren()) {
            markParents(key, value, child);
        }

        if (node.getMetaDataNames().contains(key)) {
            Object val = node.getMetaData(key);
            if (val.equals(value)) {
                // mark this node
                node.setMetaData("mark",true);
                Node n = node.getParent();
                while (n != null) {
                    n.setMetaData("mark",true);
                    n = n.getParent();
                }
            }
        }
    }

    public static void main(String[] args) throws IOException {

        FileWriter fileWriter = new FileWriter("FullTreeSampleTree.txt");
        PrintWriter pw = new PrintWriter(fileWriter);

        FBDTreeSimulator simulator = new FBDTreeSimulator();

        double lambda = 1;
        double mu = 0.5;
        double psi = 1.5;
        double x0 = 8;
        double rho = 0.5;


        int leaves = 0;
        Node myTree = null;

        while (leaves < 10 || leaves > 500) {
            myTree = simulator.simulate(lambda, mu, psi, x0, rho);
            leaves = myTree.getAllLeafNodes().size();
        }

        // printing full tree
        pw.println(myTree.toNewick());

        SampleTree sampleTree = new SampleTree();
        sampleTree.sampleTree(myTree);

        // printing sample tree
        pw.println(myTree.toNewick());

        pw.flush();
        pw.close();
    }

}
