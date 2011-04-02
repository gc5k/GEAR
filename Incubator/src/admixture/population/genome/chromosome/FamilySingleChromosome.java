package admixture.population.genome.chromosome;

public class FamilySingleChromosome {
	private int[][][] p_g; // parental chromosomes;
	// p_g[0] for dad, p_g[1] for mom, [2][haploid=2][num of locus per chr]
	private int[][][] o_g; // offspring chromosomes; [kid number][haploid=2][num of locus per chr]

	private double[][] diploid_imprint_p; //[2][haploid=2][indicator variable k indicates which population the parental allele if from]
	private double[][] diploid_imprint_o; //[Number kids][haploid=2][indicator variable k indicates which population the parental allele is from]

	private boolean disease_linked;
	private int chrID;

	public FamilySingleChromosome(int ci, int[][][] p, int[][][] o, boolean d, int N_anc) {
		chrID = ci;
		p_g = p;
		o_g = o;
		disease_linked = d;

		diploid_imprint_p = new double[2][N_anc];
		for(int i = 0; i < p_g.length; i++) {
			for(int j = 0; j < p_g[i][0].length; j++) {
				int org1 = (int) (p_g[i][0][j]/2);
				int org2 = (int) (p_g[i][1][j]/2);
				p_g[i][0][j] = p_g[i][0][j]%2;
				p_g[i][1][j] = p_g[i][1][j]%2;
				p_g[i][0][j]++;
				p_g[i][1][j]++;
				diploid_imprint_p[i][org1] += 0.5;
				diploid_imprint_p[i][org2] += 0.5;
			}
			for(int j = 0; j < diploid_imprint_p[i].length; j++) {
				diploid_imprint_p[i][j] /= p_g[0][0].length;
			}
		}

		diploid_imprint_o = new double[o_g.length][N_anc];
		for(int i = 0; i < o_g.length; i++) {
			for(int j = 0; j < o_g[i][0].length; j++) {
				int org1 = (int) (o_g[i][0][j]/2);
				int org2 = (int) (o_g[i][1][j]/2);
				o_g[i][0][j] = o_g[i][0][j]%2;
				o_g[i][1][j] = o_g[i][1][j]%2;
				o_g[i][0][j]++;
				o_g[i][1][j]++;
				diploid_imprint_o[i][org1] += 0.5;
				diploid_imprint_o[i][org2] += 0.5;
			}
			for(int j = 0; j < diploid_imprint_o[i].length; j++) {
				diploid_imprint_o[i][j] /= o_g[0][0].length;
			}
		}
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

	public int getChrID() {
		return chrID;
	}

	public double[][] getParentChromosomeAncestry() {
		return diploid_imprint_p;
	}

	public double[][] getOffspringChromosomeAncestry() {
		return diploid_imprint_o;
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
