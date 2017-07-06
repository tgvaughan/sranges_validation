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
 * Draft implementation of SR tree density described by eq 8 in Stadler et al. (2017).
 */
public class SRTreeDensity extends TreeDistribution {

    public Input <RealParameter> lambdaInput = new Input<>("lambda",
            "Speciation rate.",
            Input.Validate.REQUIRED);

    public Input <RealParameter> muInput = new Input<>("mu",
            "Extinction rate.",
            Input.Validate.REQUIRED);

    public Input <RealParameter> psiInput = new Input<>("psi",
            "Fossilization rate.",
            Input.Validate.REQUIRED);

    public Input <RealParameter> rhoInput = new Input<>("rho",
            "Present-day sampling probability.",
            Input.Validate.REQUIRED);

    public Input <RealParameter> x0Input = new Input<>("x0",
            "Time of origin of speciation process.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> sRangesInput = new Input<>("sranges",
            "Parameter representing stratigraphic ranges.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> unobsSpecTimesInput = new Input<>("unobsSpecTimes",
            "Unobserved speciation event times.",
            Input.Validate.REQUIRED);

    protected RealParameter sRanges, unobsSpecTimes;
    protected RealParameter lambda, mu, psi, rho, x0;

    protected TreeInterface tree;

    protected static final double EPSILON = 1e-10;

    public SRTreeDensity() {
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
        sRanges = sRangesInput.get();
        unobsSpecTimes = unobsSpecTimesInput.get();

    }

    protected double log_q(double t, double c1, double c2) {
        return Math.log(4) - c1*t - 2*Math.log(Math.exp(-c1*t)*(1-c2) + (1+c2));
    }

    protected double log_p(double t, double c1, double c2) {
        return Math.log(1.0 + 0.5*(-(lambda.getValue()-mu.getValue()-psi.getValue()) +
                c1*(Math.exp(-c1*t)*(1-c2)-(1+c2))/(Math.exp(-c1*t)*(1-c2)+(1+c2)))/lambda.getValue());
    }

    protected double log_qtilde_asym(double t, double c1, double c2) {
        return 0.5*(-t*(lambda.getValue()+mu.getValue()+psi.getValue()) + log_q(t, c1, c2));
    }

    protected double log_qhat_asym_ns(double startTime, double endTime, double c1, double c2) {
        return log_q(startTime, c1, c2) - log_q(endTime, c1, c2);
    }

    protected double log_qhat_asym_s(double startTime, double endTime, double c1, double c2) {
        return log_qtilde_asym(startTime, c1, c2) - log_qtilde_asym(endTime, c1, c2);
    }

    protected double log_oneMinusP0Hat(double t, double c1, double c2) {
        return Math.log(rho.getValue()*(lambda.getValue()-mu.getValue())/(lambda.getValue()*rho.getValue() + (lambda.getValue()*(1-rho.getValue()) - mu.getValue())* Math.exp((mu.getValue()-lambda.getValue()) * t))) ;
    }

    protected int enclosingStratigraphicRange(Node node) {
        Node ld = getLeftDescendantLeaf(node);
        if (node.getHeight() < (ld.getHeight() + sRanges.getArrayValue(ld.getNr())))
            return ld.getNr();
        else
            return -1;
    }

    protected Node getLeftDescendantLeaf(Node node) {
        return node.isLeaf() ? node : getLeftDescendantLeaf(node.getLeft());
    }

    @Override
    public double calculateLogP() {
        logP = 0.0;

        double c1 = Math.abs(Math.sqrt(Math.pow(lambda.getValue()-mu.getValue()-psi.getValue(),2)+4*lambda.getValue()*psi.getValue()));
        double c2 = -(lambda.getValue()-mu.getValue()-2*lambda.getValue()*rho.getValue()-psi.getValue())/c1;

        int n = tree.getLeafNodeCount();

     //   System.out.println("logP="+logP + " c1=" + c1 + " c2=" + c2);

        for (Node node : tree.getNodesAsArray()) {
            if (node.isFake())
                continue;

        //    System.out.println("node " + node.getNr());
            int enclosingRange = enclosingStratigraphicRange(node);

            double startTime, endTime;
            if (enclosingRange>=0) {

                double rangeEndTime = sRanges.getArrayValue(enclosingRange) + tree.getNode(enclosingRange).getHeight();

                // Stratigraphic range branch
                endTime = node.getHeight();
                startTime = node.isRoot() ? rangeEndTime : Math.min(rangeEndTime, node.getParent().getHeight());


                logP += log_qhat_asym_s(startTime, endTime, c1, c2);
                //System.out.println("  stratigraphic range branch("+startTime + ", " + endTime + ")");
                //System.out.println("  logP="+logP);

                // sampling contribution for oldest fossil in this stratigraphic range
                logP += Math.log(psi.getValue());

                // Non-stratigraphic range  branch

                endTime = startTime;
                startTime = node.isRoot() ? x0.getValue() : node.getParent().getHeight();

                logP += log_qhat_asym_ns(startTime, endTime, c1, c2);
                //System.out.println("  non-stratigraphic range branch("+startTime + ", " + endTime + ")");
                //System.out.println("  logP="+logP);


            } else {


                // Non-stratigraphic  branch

                endTime = node.getHeight();
                startTime = node.isRoot() ? x0.getValue() : node.getParent().getHeight();

                logP += log_qhat_asym_ns(startTime, endTime, c1, c2);
              //  System.out.println("  non-stratigraphic range branch("+startTime + ", " + endTime + ")");
              //  System.out.println("  logP="+logP);

            }

            if (node.isLeaf()) {
                if (Math.abs(node.getHeight())<EPSILON)
                    logP += Math.log(rho.getValue());
                else {
                    logP += Math.log(psi.getValue());

                    if (!node.isDirectAncestor()) {

                        logP += log_p(node.getHeight(), c1, c2);
                    //    System.out.println("  prob of survival to leaf " + node.getNr());
                    //    System.out.println("  logP="+logP);

                    } else {
                        //System.out.println("SHOULDN'T GET HERE in simple example :)");

                        // Check for descendant SR

                        Node leftDescendant = getLeftDescendantLeaf(node.getParent());
                        if (leftDescendant.getHeight() + sRanges.getValue(leftDescendant.getNr())< node.getHeight()) {
                            // We are ancestral to a distinct species: take the unobserved speciation event into account
                        }
                    }
                }
            } else {
                logP += Math.log(lambda.getValue());
            }
        }

        for (int i = 0; i < unobsSpecTimes.getDimension(); i++) {
            logP += Math.log(lambda.getValue()) + log_p(unobsSpecTimes.getArrayValue(i), c1, c2);
        }

        logP -= log_oneMinusP0Hat(x0.getValue(),c1,c2);

        return logP;
    }

       public void probSpecialCase(double y1, double y2, double o1, double o2, double x1, double x0) {

            int k = 3;

            double c1;
            double c2;
          //  double logPBis;

            double lambdaBis = lambda.getValue();
            double rhoBis = rho.getValue();
            double muBis = mu.getValue();


            System.out.println("psi\tlogP");
            for (double psiBis = 0.1; psiBis < 8.01; psiBis += 0.05) {

                psi.setValue(psiBis);

                c1 = Math.abs(Math.sqrt(Math.pow(lambdaBis - muBis - psiBis, 2) + 4 * lambdaBis * psiBis));
                c2 = -(lambdaBis - muBis - 2 * lambdaBis * rhoBis - psiBis) / c1;

                // 3 sampling events
                double logP = k * Math.log(psiBis);

                // 1 extant species
                logP += Math.log(rhoBis);

                // 1 speciation event
                logP += Math.log(lambdaBis);

                // probability of survival to present from time 3
                logP -= log_oneMinusP0Hat(3, c1, c2);

                logP += log_p(y2, c1, c2);

                logP += log_qtilde_asym(1, c1, c2) - log_qtilde_asym(0, c1, c2);


                logP += log_qtilde_asym(1.5, c1, c2) - log_qtilde_asym(0.5, c1, c2);

                logP += log_q(2,c1,c2)  + log_q(3,c1,c2) - log_q(1,c1,c2) - log_q(1.5,c1,c2);
                System.out.println(psiBis + "\t" + logP);
            }
        }

        public static void main(String[] args) {
            double y1 = 0;
            double y2 = 0.5;
            double o1 = 1;
            double o2 = 1.5;
            double x1 = 2;
            double x0 = 3;

            double lambdaBis = 1.0;
            double rhoBis = 0.1;
            double muBis = 0.1;
            double psiBis = 0.1;

            RealParameter psiParameter = new RealParameter(new Double[] {psiBis});

            List<StateNode> stateNodeList = new ArrayList<StateNode>();
            stateNodeList.add(psiParameter);

            State state = new State();
            state.stateNodeInput.setValue(stateNodeList, state);
            state.initialise();

            SRTreeDensity density = new SRTreeDensity();
            density.lambdaInput.setValue(new RealParameter(new Double[] {lambdaBis}), density);
            density.rhoInput.setValue(new RealParameter(new Double[] {rhoBis}), density);
            density.muInput.setValue(new RealParameter(new Double[] {muBis}), density);
            density.psiInput.setValue(psiParameter, density);
            density.initAndValidate();

            density.probSpecialCase(y1, y2, o1, o2, x1, x0);
        }
}
