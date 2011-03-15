package admixture.chromosome;

public class FamilySingleChromosome {
	private int[][][] p_g; // parental chromosomes;
	// p_g[0] for dad, p_g[1] for mom, [2][haploid=2][num of locus per chr]
	private int[][][] o_g; // offspring chromosomes; [kid number][haploid=2][num of locus per chr]

	private double[][][] ancestry_haploid_p; //[2][haploid=2][ancestry=num of ancetral populations]
	private double[][][] ancestry_haploid_o; //[kid number][haploid=2][ancestry=num of ancetral populations]

	private double[][] ancestry_diploid_p; //[2][ancestry=num of ancetral populations]
	private double[][] ancestry_diploid_o; //[kin number][ancestry=num of ancetral populations]

	private boolean disease_linked;
	private int chrID;
	public FamilySingleChromosome(int ci, int[][][] p, int[][][] o, boolean d) {
		chrID = ci;
		p_g = p;
		o_g = o;
		disease_linked = d;
	}

	public void AddFatherChr(int[][] g) {
		p_g[0] = g;
	}

	public void AddMotherChr(int[][] g) {
		p_g[1] = g;
	}

	public boolean isDiseaseLinked() {
		return disease_linked;
	}

	public void AscertainParentSingleChromosomeAncestry(double[][][] post_snp_ancestry) {
		ancestry_haploid_p = new double[2][2][];
		for (int i = 0; i < p_g.length; i++) {
			for (int j = 0; j < p_g[i].length; j++) {
				ancestry_haploid_p[i][j] = AscertainHaploidAncestry(p_g[i][j], post_snp_ancestry);
			}
		}
		
		ancestry_diploid_p = new double[2][ancestry_haploid_p[0][0].length];
		for (int i = 0; i < ancestry_haploid_p.length; i++) {
			for (int j = 0; j < ancestry_haploid_p[i][0].length; j++) {
				ancestry_diploid_p[i][j] = (ancestry_haploid_p[i][0][j] + ancestry_haploid_p[i][1][j])/2;
			}
		}
	}

	public void AscertainOffspringSingleChromosomeAncestry(double[][][] post_snp_ancestry) {
		ancestry_haploid_o = new double[o_g.length][2][];
		for (int i = 0; i < o_g.length; i++) {
			for (int j = 0; j < o_g[i].length; j++) {
				ancestry_haploid_o[i][j] = AscertainHaploidAncestry(o_g[i][j], post_snp_ancestry);
			}
		}

		ancestry_diploid_o = new double[ancestry_haploid_o.length][ancestry_haploid_o[0][0].length];
		for (int i = 0; i < ancestry_haploid_o.length; i++) {
			for (int j = 0; j < ancestry_haploid_o[i][0].length; j++) {
				ancestry_diploid_o[i][j] = (ancestry_haploid_o[i][0][j] + ancestry_haploid_o[i][1][j])/2;
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

	public int getChrID() {
		return chrID;
	}
	public double[][] getParentChromosomeAncestry() {
		return ancestry_diploid_p;
	}

	public double[][] getOffspringChromosomeAncestry() {
		return ancestry_diploid_o;
	}
	
	public int[] ParentGenotype(int pi, int loci ) {
		int[] d = new int[2];
		d[0] = p_g[pi][0][loci];
		d[1] = p_g[pi][1][loci];
		return d;
	}

	public int[] OffspringGenotype(int oi, int loci ) {
		int[] d = new int[2];
		d[0] = o_g[oi][0][loci];
		d[1] = o_g[oi][1][loci];
		return d;
	}

	public int[][][] getParentChromosome() {
		return p_g;
	}

	public int[][] getFatherChromosome() {
		return p_g[0];
	}

	public String getStringParentChromosome(int idx) {
		StringBuffer sb = new StringBuffer(" ");
		for(int i = 0; i < p_g[idx][0].length; i++) {
			if(p_g[idx][0][i] > p_g[idx][1][i]) {
				sb.append(p_g[idx][0][i] + " " + p_g[idx][1][i] + " ");
			} else {
				sb.append(p_g[idx][1][i] + " " + p_g[idx][0][i] + " ");
			}
		}
		return sb.toString();
	}

	public String getGenotypeStringParentChromosome(int idx) {
		StringBuffer sb = new StringBuffer(" ");
		for(int i = 0; i < p_g[idx][0].length; i++) {
			sb.append(p_g[idx][0][i] + p_g[idx][1][i] + " ");
		}
		return sb.toString();
	}

	public String getStringOffspringChromosome(int idx) {
		StringBuffer sb = new StringBuffer(" ");
		for(int i = 0; i < o_g[idx][0].length; i++) {
			if(o_g[idx][0][i] > o_g[idx][1][i]) {
				sb.append(o_g[idx][0][i] + " " + o_g[idx][1][i] + " ");
			} else {
				sb.append(o_g[idx][1][i] + " " + o_g[idx][0][i] + " ");
			}
		}
		return sb.toString();
	}

	public String getGenotypeStringOffspringChromosome(int idx) {
		StringBuffer sb = new StringBuffer(" ");
		for(int i = 0; i < o_g[idx][0].length; i++) {
			sb.append(o_g[idx][0][i] + o_g[idx][1][i] + " ");
		}
		return sb.toString();
	}

	public int[][] getMotherChromosome() {
		return p_g[1];
	}

	public int[][][] getOffspingChromosome() {
		return o_g;
	}
}
