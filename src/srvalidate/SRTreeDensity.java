package srvalidate;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeDistribution;
import beast.evolution.tree.TreeInterface;

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

        for (Node node : tree.getNodesAsArray()) {
            if (node.isFake())
                continue;

            int enclosingRange = enclosingStratigraphicRange(node);

            double startTime, endTime;
            if (enclosingRange>0) {

                double rangeEndTime = sRanges.getArrayValue(enclosingRange) + tree.getNode(enclosingRange).getHeight();

                // Stratigraphic range branch
                startTime = node.getHeight();
                endTime = node.isRoot() ? rangeEndTime : Math.min(rangeEndTime, node.getParent().getHeight());

                logP += log_qhat_asym_s(startTime, endTime, c1, c2);

                // Non-stratigraphic range  branch

                startTime = endTime;
                endTime = node.isRoot() ? x0.getValue() : node.getParent().getHeight();

                logP += log_qhat_asym_ns(startTime, endTime, c1, c2);

            } else {

                // Non-stratigraphic  branch

                startTime = node.getHeight();
                endTime = node.isRoot() ? x0.getValue() : node.getParent().getHeight();

                logP += log_qhat_asym_ns(startTime, endTime, c1, c2);
            }

            if (node.isLeaf()) {
                if (Math.abs(node.getHeight())<EPSILON)
                    logP += Math.log(rho.getValue());
                else {
                    logP += Math.log(psi.getValue());

                    if (!node.isDirectAncestor()) {
                        logP += log_p(node.getHeight(), c1, c2);

                    } else {
                        // Check for descendant SR

                        Node leftDescendant = getLeftDescendantLeaf(node.getParent());
                        if (leftDescendant.getHeight() + sRanges.getValue(leftDescendant.getNr())< node.getHeight()) {
                            // We are ancestral to a distinct species: take the unobserved speciation event into account

                            logP += Math.log(psi.getValue())
                                    + log_p(unobsSpecTimes.getArrayValue(leftDescendant.getNr()), c1, c2);
                        }
                    }
                }
            } else {
                logP += Math.log(lambda.getValue());
            }
        }

        logP -= Math.log(1.0 - log_p(x0.getArrayValue(), c1, c2));

        return logP;
    }
}
