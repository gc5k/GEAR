package admixture;

import java.util.Arrays;
import jsc.distributions.Binomial;

public class DNAStirrer {

	private int N_snp; //number of snps
	private int N_pop; //number of populations
	private int N_allele; //n_allele = 2 for a biallelic locus
	private boolean DEBUG = true;
	private int model;
	private int N; //population size
	private int round;
	private double[][] ancestral_snp_freq;
	private double[] curr_snp_freq;
	private double[] curr_snp_var;
	private String[] snp_Name;

	private double[][] pop_dna_org;
	private double[] curr_dna_org;

	private double[] pop_proportion;

	private double[][][] post_snp_prob;
	private boolean genetic_drift;
	private Binomial rnd;
	private long seed = 100;

	public DNAStirrer(int m, int s, int r, boolean gd) {
		N = s;
		model = m;
		round = r;
		genetic_drift = gd;

		N_snp = 10;
		N_pop = 1;
		N_allele = 2;
		initial_test();
	}

	private void initial_test() {
		ancestral_snp_freq = new double[N_pop][N_snp];
		for (int i = 0; i < N_pop; i++) {
			Arrays.fill(ancestral_snp_freq[i], 0.15 * (i + 1));
		}

		curr_snp_freq = new double[N_snp];
		curr_snp_var = new double[N_snp];
		curr_dna_org = new double[N_pop];

		pop_dna_org = new double[N_pop][N_pop];
		for (int i = 0; i < pop_dna_org.length; i++) {
			pop_dna_org[i][i] = 1;
		}

		pop_proportion = new double[N_pop];
		pop_proportion[0] = 1;
//		pop_proportion[1] = 0.2;
//		pop_proportion[2] = 0.025;

		post_snp_prob = new double[N][N_allele][N_pop];
	}

	public void DNAStir() {
		for (int i = 0; i < round; i++) {
			visaApplication();// migrating law may change
			double[][] M = matingMatrix();

			curr_snp_freq = new double[curr_snp_freq.length];
			curr_dna_org = new double[curr_dna_org.length];
			curr_snp_var = new double[curr_snp_var.length];

			//update allele frequency
			for (int j = 0; j < curr_snp_freq.length; j++) {
				for (int k = 0; k < M.length; k++) {
					for (int l = 0; l < M[k].length; l++) {
						curr_snp_freq[j] += M[k][l]
								* (ancestral_snp_freq[k][j] * ancestral_snp_freq[l][j] 
								+ 0.5 * ancestral_snp_freq[k][j] * (1 - ancestral_snp_freq[l][j])
								+ 0.5 * (1 - ancestral_snp_freq[k][j]) * ancestral_snp_freq[l][j]);
					}
				}
			}

			//update genomic composition
			for (int j = 0; j < curr_dna_org.length; j++) {
				for (int k = 0; k < M.length; k++) {
					for (int l = 0; l < M[k].length; l++) {
						curr_dna_org[j] += M[k][l]
								* (pop_dna_org[k][j] + pop_dna_org[l][j]) / 2;
					}
				}
			}

			//update variance for alleles
			for (int j = 0; j < curr_snp_var.length; j++) {
				for (int k = 0; k < ancestral_snp_freq.length; k++) {
					curr_snp_var[j] += curr_dna_org[k] * ancestral_snp_freq[k][j] * (1 -  ancestral_snp_freq[k][j]) / (2 * N);
				}
				curr_snp_var[j] = Math.sqrt(curr_snp_var[j]);
			}

			if(i > 0 && genetic_drift) geneticDrift();
			
			if (DEBUG) {
				System.out.println("round " + i + " current snp pool:");
				for (int j = 0; j < curr_snp_freq.length; j++) {
					System.out.print(curr_snp_freq[j] + " ");
				}
				System.out.println();

				for (int j = 0; j < curr_snp_var.length; j++) {
					System.out.print(curr_snp_var[j] + " ");
				}
				System.out.println();

				System.out.println("round " + i + " current dna composition:");
				for (int j = 0; j < curr_dna_org.length; j++) {
					System.out.print(curr_dna_org[j] + " ");
				}
				System.out.println();
			}
			System.arraycopy(curr_snp_freq, 0, ancestral_snp_freq[0], 0, curr_snp_freq.length);
			System.arraycopy(curr_dna_org, 0, pop_dna_org[0], 0, curr_dna_org.length);			
		}
	}

	private void visaApplication() {

		//update pop_proportion here

		// the criteria for application
		// 1: constant migrating rate
		for(int i = 0; i < pop_proportion.length; i++) {
			pop_proportion[i] = pop_proportion[i];
		}
		// 2: policy varies
		// for marry in model, the 2p1-1>0
	}

	private void geneticDrift() {
		rnd = new Binomial(2*N, 0.5);
		rnd.setSeed(seed);
		for(int i = 0; i < curr_snp_freq.length; i++) {
			rnd.setP(curr_snp_freq[i]);
			curr_snp_freq[i] = rnd.random()/(2*N);
		}
	}

	private double[][] matingMatrix() {
		double[][] M = null;
		if (model == AdmixtureConstant.Marry_In) {
			M = new double[pop_proportion.length][pop_proportion.length];
			M[0][0] = 1 - 2 * (1 - pop_proportion[0]);
			for (int i = 1; i < M[0].length; i++) {
				M[0][i] = pop_proportion[i];
				M[i][0] = pop_proportion[i];
			}
		} else if (model == AdmixtureConstant.Mirgrate_In) {
			M = Tools.Generate_Matrix(pop_proportion);
		}
		return M;
	}

	public double[][] get_Frequency() {
		return ancestral_snp_freq;
	}

	public static void main(String[] args) {

		DNAStirrer ds = new DNAStirrer(2, 10000, 2, true);
		ds.DNAStir();

		double[][] sf = ds.get_Frequency();

		System.out.println("========");
		for (int i = 0; i < sf.length; i++) {
			for (int j = 0; j < sf.length; j++) {
				System.out.print(sf[i][j] + " ");
			}
			System.out.println();
		}
	}

}
