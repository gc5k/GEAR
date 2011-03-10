package admixture;

import java.util.Arrays;

import arsenal.Tools;
import jsc.distributions.Binomial;

/**
*
* @author Guo-Bo Chen, chenguobo@gmail.com
*/

public class DNAStirrer {

	private int N_snp; //number of snps
	private int N_pop; //number of populations
	private int N_allele; //n_allele = 2 for a biallelic locus
	private boolean DEBUG = false;
	private int model;
	private int N; //population size
	private int round;
	private double[][] ancestral_snp_freq; //this one records the frequencies for the ancestral populations
											// it should not be changed
	private double[][] src_snp_freq; //this one records the frequencies of the populations in mating;
									//the first row is the updated frequencies of the admixed population; 
	private double[] curr_snp_freq; //the mixed frequency of the current population under random mating
	private double[] curr_snp_var; //the mixed variance of the current population under random mating
	private String[] snp_Name;

	private double[][] pop_dna_org;
	private double[] curr_dna_org;

	private double[] pop_weight;

	private double[][][] post_snp_prob;//post_snp_prob[locus][allele][ancestry]
	private boolean genetic_drift;
	private Binomial rnd;
	private long seed = 100;

	public DNAStirrer(int m, int s, int ns, boolean gd) {//for test

		model = m;
		N = s;
		genetic_drift = gd;

		N_snp = ns;
		N_pop = 2;
		N_allele = 2;
		initial_test();
	}

	private void initial_test() {
		src_snp_freq = new double[N_pop][N_snp];
		for (int i = 0; i < N_pop; i++) {
			Arrays.fill(src_snp_freq[i], 0.15 * (i + 1));
		}

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

		pop_weight = new double[N_pop];
		pop_weight[0] = 0.8;
		pop_weight[1] = 0.2;
//		pop_weight[2] = 0.025;

		post_snp_prob = new double[N_snp][N_allele][N_pop];
	}

	public void DNAStir(int r) {
		round = r;
		for (int i = 0; i < round; i++) {
			double[][] M = visaOffice();

			curr_snp_freq = new double[curr_snp_freq.length];
			curr_dna_org = new double[curr_dna_org.length];
			curr_snp_var = new double[curr_snp_var.length];

			//update genomic composition
			for (int j = 0; j < curr_dna_org.length; j++) {
				for (int k = 0; k < M.length; k++) {
					for (int l = 0; l < M[k].length; l++) {
						curr_dna_org[j] += M[k][l]
								* (pop_dna_org[k][j] + pop_dna_org[l][j]) / 2;
					}
				}
			}

			//update allele frequency
			for (int j = 0; j < curr_snp_freq.length; j++) {
				for (int k = 0; k < M.length; k++) {
					for (int l = 0; l < M[k].length; l++) {
						curr_snp_freq[j] += M[k][l]
								* (src_snp_freq[k][j] * src_snp_freq[l][j] 
								+ 0.5 * src_snp_freq[k][j] * (1 - src_snp_freq[l][j])
								+ 0.5 * (1 - src_snp_freq[k][j]) * src_snp_freq[l][j]);
					}
				}
			}

			//update variance for alleles
			for (int j = 0; j < curr_snp_var.length; j++) {
				for (int k = 0; k < src_snp_freq.length; k++) {
					curr_snp_var[j] += curr_dna_org[k] * src_snp_freq[k][j] * (1 -  src_snp_freq[k][j]) / (2 * N);
				}
				curr_snp_var[j] = Math.sqrt(curr_snp_var[j]);
			}

			for (int j = 0; j < post_snp_prob.length; j++) {
				double c0 = 0;
				double c1 = 0;
				for (int k = 0; k < curr_dna_org.length; k++) {
					c0 += ancestral_snp_freq[k][j] * curr_dna_org[k];
					c1 += (1-ancestral_snp_freq[k][j]) * curr_dna_org[k];
				}
				for (int k = 0; k < curr_dna_org.length; k++) {
					post_snp_prob[j][0][k] = ancestral_snp_freq[k][j] * curr_dna_org[k]/c0;
					post_snp_prob[j][1][k] = (1 - ancestral_snp_freq[k][j]) * curr_dna_org[k]/c1;
				}
			}
			//update posterior probability that an allele is contributed by a certain population;
			
			if(i > 0 && genetic_drift) geneticDrift();
			
			if (DEBUG) {
				System.out.println("====round " + i + "\n current snp pool:");
				for (int j = 0; j < curr_snp_freq.length; j++) {
					System.out.print(curr_snp_freq[j] + " ");
				}
				System.out.println();

				for (int j = 0; j < curr_snp_var.length; j++) {
					System.out.print(curr_snp_var[j] + " ");
				}
				System.out.println();

				System.out.println(" current dna composition:");
				for (int j = 0; j < curr_dna_org.length; j++) {
					System.out.print(curr_dna_org[j] + " ");
				}
				System.out.println();
				
				System.out.println(" posterior snp probability:");
				for (int j = 0; j < post_snp_prob.length; j++) {
					System.out.println("snp:" + j);
					for (int k = 0; k < post_snp_prob[j].length; k++) {
						for (int l = 0; l < post_snp_prob[j][k].length; l++) {
							System.out.print(post_snp_prob[j][k][l] + " ");
						}
						System.out.println();
					}
				}
				System.out.println();
			}
			System.arraycopy(curr_snp_freq, 0, src_snp_freq[0], 0, curr_snp_freq.length);
			System.arraycopy(curr_dna_org, 0, pop_dna_org[0], 0, curr_dna_org.length);			
		}
	}

	private void geneticDrift() {
		rnd = new Binomial(2*N, 0.5);
		rnd.setSeed(seed);
		for(int i = 0; i < curr_snp_freq.length; i++) {
			rnd.setP(curr_snp_freq[i]);
			curr_snp_freq[i] = rnd.random()/(2*N);
		}
	}

	private double[][] visaOffice() {
		//update pop_weight here

		// the criteria for application
		// 1: constant migrating rate
		for(int i = 0; i < pop_weight.length; i++) {
			pop_weight[i] = pop_weight[i];
		}
		// 2: policy varies
		// for marry in model, the 2p1-1>0
		double[][] M = null;
		if (model == AdmixtureConstant.Marry_In) {
			M = new double[pop_weight.length][pop_weight.length];
			M[0][0] = 1 - 2 * (1 - pop_weight[0]);
			for (int i = 1; i < M[0].length; i++) {
				M[0][i] = pop_weight[i];
				M[i][0] = pop_weight[i];
			}
		} else if (model == AdmixtureConstant.Mirgrate_In) {
			M = Tools.Generate_Matrix(pop_weight);
		}
		return M;
	}

	public double[][] MultiPopSNPFrequency() {
		return src_snp_freq;
	}

	public double[] SNPPanel() {
		return curr_snp_freq;
	}

	public double[][][] PostSNPAncestralProb() {
		return post_snp_prob;
	}
	
	public int NumberOfSNP() {
		return N_snp;
	}

	public static void main(String[] args) {

		DNAStirrer ds = new DNAStirrer(2, 10000, 2, AdmixtureConstant.With_Genetic_Drift);
		ds.DNAStir(10);

		double[][] sf = ds.MultiPopSNPFrequency();

		System.out.println("========");
		for (int i = 0; i < sf.length; i++) {
			for (int j = 0; j < sf[i].length; j++) {
				System.out.print(sf[i][j] + " ");
			}
			System.out.println();
		}
	}

}
