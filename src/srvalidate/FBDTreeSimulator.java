package srvalidate;

import beast.evolution.tree.Node;
import beast.util.Randomizer;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

/**
 * A class that simulates full trees under the fossilized birth-death process with budding speciation.
 *
 * @author Alexei Drummond
 */
public class FBDTreeSimulator {

    /**
     * Simulate a full phylogeny under the fossilized birth-death process with budding speciation.
     *
     * @param lambda speciation rate
     * @param mu     extinction rate
     * @param psi    fossil sampling rate
     * @param x0     time of the origin
     * @param rho    probability of sampling an extant species at the present
     *               returned sample phylogeny.
     * @return a Node representing the root of the full phylogeny.
     */
    public Node simulate(double lambda, double mu, double psi, double x0, double rho) {


        Set<String> taxa = new TreeSet<>();

        double T = 0;

        Node start = createNode(taxa);
        start.setHeight(x0);

        List<Node> activeLineages = new ArrayList<>();

        activeLineages.add(start);

        int psiSamples = 0;
        int rhoSamples = 0;

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
                    doExtinction(activeLineages, T, x0);

                } else {
                    // do sample
                    Node sampledNode = activeLineages.get(Randomizer.nextInt(activeLineages.size()));
                    doFossilSample(sampledNode, T, x0, taxa, activeLineages);
                    psiSamples += 1;
                }
            }
        }

        for (Node activeLineage : activeLineages) {
            if (Randomizer.nextDouble() < rho) {
                doRho(activeLineage);
                rhoSamples += 1;
            } else {
                doUnsample(activeLineage);
            }
        }

        processMetaData(start);

        return start;
    }

    private void doUnsample(Node node) {
        node.setHeight(0);
        node.setMetaData("rho", false);
    }

    private void doRho(Node node) {
        node.setHeight(0);
        node.setMetaData("rho", true);
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
        parent.setMetaData("reaction", "lambda");

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
    private void doExtinction(List<Node> nodes, double time, double x0) {

        Node deadNode = nodes.get(Randomizer.nextInt(nodes.size()));
        deadNode.setHeight(x0 - time);
        deadNode.setMetaData("reaction", "mu");
        nodes.remove(deadNode);
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
    private void doFossilSample(Node nodeToSample, double time, double x0, Set<String> taxa, List<Node> nodes) {

        Node continuingLineage = createNode(taxa);

        nodeToSample.setHeight(x0 - time);
        nodeToSample.setMetaData("reaction", "psi");

        nodeToSample.addChild(continuingLineage);
        continuingLineage.setHeight(x0 - time);

        nodes.add(continuingLineage);
        nodes.remove(nodeToSample);
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

        FileWriter fileWriter = new FileWriter("FBDTreeSimulator.txt");
        PrintWriter pw = new PrintWriter(fileWriter);

        FBDTreeSimulator simulator = new FBDTreeSimulator();

        double lambda = 1;
        double mu = 0.5;
        double psi = 1.5;
        double x0 = 8;
        double rho = 0.5;

        int rep = 0;

        pw.println("n\trootHeight\tn_sample\trootHeight_sample");
        while (rep < 10000) {
            Node myTree = simulator.simulate(lambda, mu, psi, x0, rho);
            int leaves = myTree.getAllLeafNodes().size();
            double height = myTree.getHeight();
            SampleTree sampleTree = new SampleTree();
            sampleTree.sampleTree(myTree);
            int sleaves = myTree.getAllLeafNodes().size();
            double sheight = myTree.getHeight();

            if (leaves>4) {
                pw.println(leaves + "\t" + myTree.getHeight() + "\t" + sleaves + "\t" + sheight);
                rep += 1;
            }
        }
        pw.flush();
        pw.close();
    }
}
