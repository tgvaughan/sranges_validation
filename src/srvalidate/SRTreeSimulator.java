package srvalidate;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeUtils;
import beast.util.Randomizer;

import java.util.*;

/**
 * Created by lauragarcia on 26/06/17.
 */
public class SRTreeSimulator {


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
                    doDeath(activeLineages, T, x0, srangeMap);

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

        // TODO traverse the whole tree, removing subtrees that have no stratigraphic intervals in them
        for (Node node : start.getAllChildNodes()) {
            for (Node childNode : node.getAllChildNodes()) {
                if (srangeMap.get(childNode) == null) {
                    node.removeChild(childNode);
                }
            }
        }

        processMetaData(start);

        System.out.println("Random tree has " + srangeMap.keySet().size() + " stratigraphic ranges." + srangeMap.values());

        return start;
    }

    private Node createNode(Set<String> taxa) {
        String newTaxon = (taxa.size()+1)+"";
        taxa.add(newTaxon);
        Node newNode = new Node(newTaxon);
        newNode.setNr(taxa.size()-1);

        return newNode;
    }

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

    private void doDeath(List<Node> nodes, double time, double x0, Map<Node,RealParameter> srangeMap) {

        Node deadNode = nodes.get(Randomizer.nextInt(nodes.size()));
        nodes.remove(deadNode);

        //if (srangeMap.get(deadNode) == null) {
        //    Node parent = deadNode.getParent();
        //    parent.removeChild(deadNode);
        //}
    }

    private void doSample(Node nodeToSample, double time, double x0, Map<Node,RealParameter> srangeMap) {

        RealParameter srange = srangeMap.get(nodeToSample);
        if (srange == null) {
            RealParameter newsrange = new RealParameter(new Double[2]);
            newsrange.setValue(0,x0-time);
            newsrange.setValue(1,x0-time);
            srangeMap.put(nodeToSample,newsrange);
            nodeToSample.setMetaData("oldestFossil",x0-time);
        } else {
            srange.setValue(1,x0-time);
        }
        nodeToSample.setHeight(x0-time);
    }

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


    public static void main(String[] args) {

        SRTreeSimulator simulator = new SRTreeSimulator();

        Node myTree = simulator.simulate(1,0.25,1.5,3,0.5);

        String treeInNewick = myTree.toNewick();
        System.out.println(treeInNewick);
    }
}
