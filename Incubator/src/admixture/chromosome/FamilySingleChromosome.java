package admixture.chromosome;

public class FamilySingleChromosome {
	private int[][][] p_g; // parental chromosomes;
	// p_g[0] for dad, p_g[1] for mom, [2][haploid=2][locus]
	private int[][][] o_g; // offspring chromosomes; [kid number][haploid=2][locus]

	private double[][][] ancestry_haploid_p; //[2][haploid=2][ancestry]
	private double[][][] ancestry_haploid_o; //[kid number][haploid=2][ancestry]

	public FamilySingleChromosome(int[][][] p, int[][][] o) {
		p_g = p;
		o_g = o;
	}

	public FamilySingleChromosome(int id, int num_kid) {
		p_g = new int[2][][];
		o_g = new int[num_kid][][];
	}

	public void AddFatherChr(int[][] g) {
		p_g[0] = g;
	}

	public void AddMotherChr(int[][] g) {
		p_g[1] = g;
	}

	public void AscertainParentSingleChromosomeAncestry(double[][][] post_snp_ancestry) {
		ancestry_haploid_p = new double[2][2][];
		for (int i = 0; i < p_g.length; i++) {
			for (int j = 0; j < p_g[i].length; j++) {
				ancestry_haploid_p[i][j] = AscertainHaploidAncestry(p_g[i][j], post_snp_ancestry);
			}
		}
	}

	public void AscertainOffspringSingleChromosomeAncestry(double[][][] post_snp_ancestry) {
		ancestry_haploid_o = new double[2][2][];
		for (int i = 0; i < o_g.length; i++) {
			for (int j = 0; j < o_g[i].length; j++) {
				ancestry_haploid_o[i][j] = AscertainHaploidAncestry(o_g[i][j], post_snp_ancestry);
			}
		}		
	}

	private double[] AscertainHaploidAncestry(int[] g, double[][][] post_snp_ancestry) {
		double[] ancestry = new double[post_snp_ancestry[0][0].length];
		for (int i = 0; i < g.length; i++) {//allele
			for (int j = 0; j < post_snp_ancestry[i][g[i]].length; j++) {//ancestry
				ancestry[j] += post_snp_ancestry[i][g[i]][j];
			}
		}
		for (int i = 0; i < ancestry.length; i++) {
			ancestry[i] /= g.length;
		}
		return ancestry;
	}
	
	public double[][] getParentSingleChromosomeAncestry() {
		double[][] a_p = new double[ancestry_haploid_p.length][ancestry_haploid_p[0][0].length];
		for (int i = 0; i < 0; i++) {
			for (int j = 0; j < a_p[i].length; i++) {
				a_p[i][j] = (ancestry_haploid_p[i][0][j] + ancestry_haploid_p[i][1][j])/2;
			}
		}
		return a_p;
	}
	
	public double[][] getOffspringSingleChromosomeAncestry() {
		double[][] a_o = new double[ancestry_haploid_p.length][ancestry_haploid_p[0][0].length];
		for (int i = 0; i < 0; i++) {
			for (int j = 0; j < a_o[i].length; i++) {
				a_o[i][j] = (ancestry_haploid_o[i][0][j] + ancestry_haploid_o[i][1][j])/2;
			}
		}
		return a_o;
	}
}
