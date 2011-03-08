package admixture;

public class SequenceGenerator {

	private double[] snp_panel;
	private double[][][] posterior_probability;

	public SequenceGenerator(double[] sp, double[][][] pp) {
		snp_panel = sp;
		posterior_probability = pp;
	}
}
