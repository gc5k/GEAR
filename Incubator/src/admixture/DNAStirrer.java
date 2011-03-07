package admixture;

import java.util.Arrays;
import org.apache.commons.lang3.ArrayUtils;
public class DNAStirrer {

	private int model;
	private int population_size;
	private int round;
	private double[][] snp_Freq;
	private double[] curr_snp_pool;
	private double[][] pop_dna_composition;
	private double[] curr_dna_composition;
	private String[] snp_Name;

	private double[] population_proportion;
	
	public DNAStirrer(int m, int s, int r) {
		model = m;
		population_size = s;
		round = r;
		snp_Freq = new double[3][10];
		for(int i = 0; i < snp_Freq.length; i++) {
			Arrays.fill(snp_Freq[i], 0.15*(i+1));
		}
		curr_snp_pool = new double[snp_Freq[0].length];
		curr_dna_composition = new double[snp_Freq.length];
		
		pop_dna_composition = new double[snp_Freq.length][snp_Freq.length];
		for(int i = 0; i < pop_dna_composition.length; i++) {
			pop_dna_composition[i][i] = 1;
		}

		population_proportion = new double[snp_Freq.length];
		population_proportion[0] = 0.7;
		population_proportion[1] = 0.2;
		population_proportion[2] = 0.1;
	}
	
	public void DNAStir() {
		for(int i = 0; i < round; i++) {
			double[][] M = matingMatrix();
			
			curr_snp_pool = new double[curr_snp_pool.length];
			curr_dna_composition = new double[curr_dna_composition.length];
			for(int j = 0; j < curr_snp_pool.length; j++) {
				for(int k = 0; k < M.length; k++) {
					for(int l = 0; l < M[k].length; l++) {
						curr_snp_pool[j] += M[k][l] * (snp_Freq[k][j] * snp_Freq[l][j] 
						                 + 0.5 * snp_Freq[k][j] * (1 - snp_Freq[l][j]) 
						                 + 0.5 * (1 - snp_Freq[k][j]) * snp_Freq[l][j] );
					}
				}
			}
			System.out.println("current snp pool:");
			for(int j = 0; j < curr_snp_pool.length; j++) {
				System.out.print(curr_snp_pool[j] + " ");
			}
			System.out.println("\n");

			for(int j = 0; j < curr_dna_composition.length; j++) {
				for(int k = 0; k < M.length; k++) {
					for(int l = 0; l < M[k].length; l++) {
						curr_dna_composition[j] += M[k][l] * (pop_dna_composition[k][j] + pop_dna_composition[l][j])/2;
					}
				}			
			}

			System.out.println("current dna composition:");
			for(int j = 0; j < curr_dna_composition.length; j++) {
				System.out.print(curr_dna_composition[j] + " ");
			}
			System.out.println("\n");
			System.out.println();
			System.arraycopy(curr_snp_pool, 0, snp_Freq[0], 0, curr_snp_pool.length);
			System.arraycopy(curr_dna_composition, 0, pop_dna_composition[0], 0, curr_dna_composition.length);
		}
	}
	
	private double[][] matingMatrix() {
		double[][] M = null;
		if(model == AdmixtureConstant.Marry_In) {
			M = new double[population_proportion.length][population_proportion.length];
			M[0][0] = 1 - 2 * (1 - population_proportion[0]);
			for (int i = 1; i < M[0].length; i++) {
				M[0][i] = population_proportion[i];
				M[i][0] = population_proportion[i];
			}
		} else if (model == AdmixtureConstant.Mirgrate_In) {
			M = Tools.Generate_Matrix(population_proportion);
		}
		return M;
	}

	public double[][] get_Frequency() {
		return snp_Freq;
	}

	public static void main(String[] args) {

		DNAStirrer ds = new DNAStirrer(1, 100, 2);
		ds.DNAStir();
		double[][] sf = ds.get_Frequency();
		for(int i = 0; i < sf.length; i++) {
			for(int j = 0; j < sf.length; j++) {
				System.out.print(sf[i][j] + " ");
			}
			System.out.println();
		}
	}

}
