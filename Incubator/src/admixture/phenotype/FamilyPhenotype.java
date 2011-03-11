package admixture.phenotype;

public class FamilyPhenotype {
	private double[][] p_phe;
	private int[] p_status;
	private double[][] o_phe;
	private int[] o_status;

	private int num_kid;
	public FamilyPhenotype(double[][] p_p, int[] p_s, double[][] o_p, int[] o_s ) {
		p_phe = p_p;
		p_status = p_s;
		o_phe = o_p;
		o_status = o_s;
	}

	public void GeneratePhenotypeNull_I() {
		
	}
}
