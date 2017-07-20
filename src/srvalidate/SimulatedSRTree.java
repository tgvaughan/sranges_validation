package srvalidate;

import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Alexei Drummond
 */
public class SimulatedSRTree extends Tree {

    public Input<RealParameter> lambdaInput = new Input<>("lambda",
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

    public Input <IntegerParameter> minSizeInput = new Input<>("minSize",
            "The minimum number of tips for the tree.",
            Input.Validate.REQUIRED);

    public void initAndValidate() {

        SRTreeSimulator simulator = new SRTreeSimulator();

        double lambda = lambdaInput.get().getValue();
        double mu = muInput.get().getValue();
        double psi = psiInput.get().getValue();
        double x0 = x0Input.get().getValue();
        double rho = rhoInput.get().getValue();

        int n = 0;

        List<Double> unobservedSpeciationAges = new ArrayList<>();
        List<Double> stratigraphicIntervals = new ArrayList<>();
        Node myTree = null;
        while (n < minSizeInput.get().getValue()) {
            unobservedSpeciationAges.clear();
            stratigraphicIntervals.clear();
            myTree = simulator.simulate(lambda, mu, psi, x0, rho, unobservedSpeciationAges, stratigraphicIntervals);
            n = myTree.getAllLeafNodes().size();
        }

        setRoot(myTree);
    }

}
