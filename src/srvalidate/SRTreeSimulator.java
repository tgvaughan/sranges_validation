package srvalidate;

import beast.evolution.tree.Node;
import beast.util.Randomizer;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by lauragarcia on 26/06/17.
 */
public class SRTreeSimulator {


    public Node simulate(double lambda, double mu, double psi, double x0) {

        double T = 0;

        Node start = new Node();
        start.setHeight(x0);

        List<Node> unsampled = new ArrayList<>();
        List<Node> sampled = new ArrayList<>();

        unsampled.add(start);

        while (T < x0) {

            int n = unsampled.size() + sampled.size();

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

                // draw random
                // calculate if unsampled parent
                if (true) { // unsampled

                    int grandparent = Randomizer.nextInt(unsampled.size());

                    Node parent = new Node();
                    parent.setHeight(x0 - T);
                    parent.setParent(unsampled.get(grandparent));
                    Node rightChild = new Node();
                    rightChild.setHeight(x0 - T);
                    parent.addChild(rightChild);

                    unsampled.remove(grandparent);
                    unsampled.add(parent);
                    unsampled.add(rightChild);

                } else { // sampled parent

                }

            } else {
                U -= totalBirthRate/totalPropensity;
                if (U < totalDeathRate/totalPropensity) {
                    // do death
                } else {
                    // do sample
                }
            }
        }

    }
}
