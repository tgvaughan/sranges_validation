package srvalidate;

import beast.evolution.tree.Node;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * Created by lauragarcia on 25/07/17.
 */
public class SRTree {

    public SRTree() {}

        /**
         * Takes a sample tree and produce the stratigraphic intervals
         * @param sampleTree
         */
        public static void srTree(Node sampleTree) {

            for (Node node : sampleTree.getChildren()) {
                boolean isAlreadySampled = false;
                if (node.getMetaData("reaction")=="psi") {
                    if (isAlreadySampled = false)
                        node.setMetaData("oldestFossil", node.getHeight());
                    else
                        node.setMetaData("youngestFossil", node.getHeight());
                    isAlreadySampled = true;


                }

            }
        }

    public static void main(String[] args) throws IOException {

        FileWriter fileWriter = new FileWriter("FullTreeSampleTreeSRTree.txt");
        PrintWriter pw = new PrintWriter(fileWriter);

        FBDTreeSimulator simulator = new FBDTreeSimulator();

        double lambda = 1;
        double mu = 0.5;
        double psi = 1.5;
        double x0 = 8;
        double rho = 0.5;


        int leaves = 0;
        Node myTree = null;

        while (leaves < 10) {
            myTree = simulator.simulate(lambda, mu, psi, x0, rho);
            leaves = myTree.getAllLeafNodes().size();
        }

        // printing full tree
        pw.println(myTree.toNewick());

        SampleTree sampleTree = new SampleTree();
        sampleTree.sampleTree(myTree);

        // printing sample tree
        pw.println(myTree.toNewick());

        SRTree srTree = new SRTree();
        SRTree.srTree(myTree);

        // printing sr tree
        pw.println(myTree.toNewick());

        pw.flush();
        pw.close();
    }
    }

