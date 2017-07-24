package srvalidate;

import beast.evolution.tree.Node;

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

        removeNodesByMetadata("reaction", "mu", fullTree);
        removeNodesByMetadata("rho", false, fullTree);
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
}
