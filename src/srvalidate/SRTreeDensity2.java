package srvalidate;

import beast.core.Input;
import beast.core.State;
import beast.core.StateNode;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeDistribution;
import beast.evolution.tree.TreeInterface;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by lauragarcia on 8/08/17.
 */


public class SRTreeDensity2 extends TreeDistribution {

    public Input<RealParameter> lambdaInput = new Input<>("lambda",
            "Speciation rate.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> muInput = new Input<>("mu",
            "Extinction rate.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> psiInput = new Input<>("psi",
            "Fossilization rate.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> rhoInput = new Input<>("rho",
            "Present-day sampling probability.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> x0Input = new Input<>("x0",
            "Time of origin of speciation process.",
            Input.Validate.REQUIRED);


    protected RealParameter lambda, mu, psi, rho, x0;

    protected TreeInterface tree;

    public SRTreeDensity2() {
        treeInput.setRule(Input.Validate.REQUIRED);
        treeIntervalsInput.setRule(Input.Validate.FORBIDDEN);
    }

    @Override
    public void initAndValidate() {
        tree = treeInput.get();

        lambda = lambdaInput.get();
        mu = muInput.get();
        psi = psiInput.get();
        rho = rhoInput.get();
        x0 = x0Input.get();

    }

    public double p(double t, double c1, double c2) {
        return 1 + ((-(lambda.getValue() - mu.getValue() - psi.getValue()) + c1 * ((Math.exp(-c1 * t) * (1 - c2) - (1 + c2)) / (Math.exp(-c1 * t) * (1 - c2) + (1 + c2)))) / (2 * lambda.getValue()));
    }

    public double q(double t, double c1, double c2) {
        return (4 * Math.exp(-c1 * t)) / (Math.exp(-c1 * t) * (1 - c2) + Math.pow(1 + c2, 2));
    }

    public double q_tilde_asym(double t, double c1, double c2) {
        return Math.sqrt(Math.exp(-t * (lambda.getValue() + mu.getValue() + psi.getValue())) * q(t, c1, c2));
    }

    public double P(double t, double c1, double c2) {
        return ((lambda.getValue() + mu.getValue() + psi.getValue() - c1) * t - 2 * Math.log(Math.exp(-c1 * t) * (1 - c2) + (1 + c2))) / (2 * lambda.getValue());
    }


    /**
     * Modifies logP according to the contribution of the set P
     *
     * @param i,j,logP,c1,c2
     */
    private void setPContribution(Node i, Node j, double logP, double c1, double c2) {

        if (j.getMetaDataNames().contains(SampleTree.S_RANGE)) {

            // if stratigraphic range i is ancestral to stratigraphic range j
            if (i.getHeight() < (double) j.getMetaData(SampleTree.S_RANGE) + j.getHeight()) {

                // if it exists a stratigraphic range ancestral to stratigraphic range j which is not also ancestral to stratigraphic range i, exit
                for (Node k : tree.getNodesAsArray()) {
                    if (k.getHeight() < (double) j.getMetaData(SampleTree.S_RANGE) + j.getHeight()) {
                        if (k.getHeight() >= (double) i.getMetaData(SampleTree.S_RANGE) + i.getHeight()) {
                            return;
                        }
                    }
                }

                // if ranges i and j are lied on a straight line
                Node node = i;
                while (node.getChildCount() > 0) {
                    if (node.equals(j)) {
                        logP += Math.log(lambda.getValue()) + Math.log(P(i.getHeight(), c1, c2) - P((double) j.getMetaData(SampleTree.S_RANGE) + j.getHeight(), c1, c2));
                        break;
                    }
                    node = node.getChild(0);
                }
            }
        }
    }

    @Override
    public double calculateLogP() {

        logP = 0.0;

        double c1 = Math.abs(Math.sqrt(Math.pow(lambda.getValue() - mu.getValue() - psi.getValue(), 2) + 4 * lambda.getValue() * psi.getValue()));
        double c2 = -(lambda.getValue() - mu.getValue() - 2 * lambda.getValue() * rho.getValue() - psi.getValue()) / c1;

        double Ls = 0.0;
        int k = 0;
        int l = 0;
        int n = tree.getLeafNodeCount();

        for (Node node : tree.getNodesAsArray()) {

            if (node.isFake())
                continue;

            //psi sampling
            if (SampleTree.isPsiSampleNode(node)) {

                // if the stratigraphic range is represented by a single fossil
                if ((double) node.getMetaData(SampleTree.S_RANGE) == 0.0) {
                    k += 1;

                    if (node.isLeaf()) {
                        logP += Math.log(p(node.getHeight(), c1, c2));
                    }
                }
                // if the stratigraphic range have different start and end times
                else {
                    double endTime = node.getHeight();
                    double startTime = node.getHeight() + (double) node.getMetaData(SampleTree.S_RANGE);

                    k += 2;
                    Ls += endTime - startTime;
                    logP += Math.log(q_tilde_asym(startTime, c1, c2)) - Math.log(q_tilde_asym(endTime, c1, c2));

                    if (node.isLeaf()) {
                        logP += Math.log(p(node.getHeight(), c1, c2));
                    }
                }

                for (Node j : tree.getNodesAsArray()) {

                    if (j == null) throw new RuntimeException("This should not be possible!");

                    setPContribution(node, j, logP, c1, c2);
                }
            }

            //rho sampling
            else {


                if (node.isLeaf() && node.getMetaData("rho").equals("true")) {
                    l += 1;
                    // if the stratigraphic range is represented by a single fossil
                    if ((double) node.getMetaData(SampleTree.S_RANGE) == 0.0) {
                        k += 1;
                    } else {
                        double endTime = node.getHeight();
                        double startTime = node.getHeight() + (double) node.getMetaData(SampleTree.S_RANGE);
                        k += 2;
                        Ls += endTime - startTime;
                        logP += Math.log(q_tilde_asym(startTime, c1, c2)) - Math.log(q_tilde_asym(endTime, c1, c2));
                    }

                    for (Node j : tree.getNodesAsArray()) {
                        setPContribution(node, j, logP, c1, c2);
                    }
                }

                // speciation
                else {
                    logP += Math.log(lambda.getValue());
                    double startTime;
                    if (node.getParent() == null) {
                        startTime = x0.getValue();
                    } else {
                        startTime = node.getParent().getHeight();
                    }

                    double endTime = node.getHeight();
                    logP += Math.log(q(startTime, c1, c2)) - Math.log(q(endTime, c1, c2));

                }
            }

        }

        logP += k * Math.log(psi.getValue());
        logP += l * Math.log(rho.getValue());
        logP += psi.getValue() * Ls;

        logP -= Math.log(1 - p(x0.getValue(), c1, c2));


        return logP;
    }

    public static void main(String[] args) {

        double lambdaBis = 1.0;
        double rhoBis = 0.5;
        double muBis = 0.5;
        double psiBis = 2.5;
        double x0Bis = 8.0;

        RealParameter lambdaParameter = new RealParameter(new Double[]{lambdaBis});

        List<StateNode> stateNodeList = new ArrayList<StateNode>();
        stateNodeList.add(lambdaParameter);

        State state = new State();
        state.stateNodeInput.setValue(stateNodeList, state);
        state.initialise();

        SRTreeDensity2 density = new SRTreeDensity2();
        density.lambdaInput.setValue(lambdaParameter, density);
        density.rhoInput.setValue(new RealParameter(new Double[]{rhoBis}), density);
        density.muInput.setValue(new RealParameter(new Double[]{muBis}), density);
        density.psiInput.setValue(new RealParameter(new Double[]{psiBis}), density);
        density.x0Input.setValue(new RealParameter(new Double[]{x0Bis}), density);
        density.initAndValidate();
    }

}