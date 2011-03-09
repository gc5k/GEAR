package admixture;

import java.util.Arrays;
import java.util.Random;
import jsc.distributions.Poisson;
import jsc.distributions.Uniform;

public class HotSpot {

	private boolean DEBUG = true;
	private double len_Morgan; // length of the chromosome measured by Morgan
	private int N_snp;
	private Random rnd = new Random();
	private double[] snp_panel;

	private double[][] rec_frac; // there two kinds of recombination fractions
	// 1: free of recombination that each element equals 0.5, rec_frac=[0.5,0.5,0.5,...]
	// it is generated in the method recombinationFree()
	// 2: cross over between snps: rec_frac looks like=[0, 0, 0, 1, 0, 0, 1, ...]
	// it is generated in the method recombination()

	private int[][] rec_hotspot;

	public HotSpot(double[] sp) {
		snp_panel = sp;

		N_snp = snp_panel.length;

		len_Morgan = 1;
		rec_frac = new double[2][N_snp];
		rec_hotspot = new int[2][];
	}

	public void GenerateChromosome(boolean rf) {

		for(int i = 0; i < 10; i++) {
			if (rf) {
				recombinationFree();
			} else {
				recombination();
			}
		}
	}

	private void recombinationFree() {

		for (int i = 0; i < rec_frac.length; i++) {
			rec_hotspot[i] = new int[1];
			rec_hotspot[i][0] = 0;
			Arrays.fill(rec_frac[i], 0.5);
		}
		if (DEBUG) print();
	}

	private void recombination() {

		Poisson p = new Poisson(len_Morgan);
		Uniform u = new Uniform(1, N_snp - 1);

		for (int i = 0; i < rec_hotspot.length; i++) {
			if (rec_hotspot[i] != null) {
				Tools.Fill_Matrix(rec_frac, i, rec_hotspot[i], 0);
			}

			int rec_event = 0;
			do {
				rec_event = (int) p.random();
			} while (rec_event >= N_snp - 2);
			rec_hotspot[i] = new int[rec_event + 1];

			int c = 1;
			while (c < rec_hotspot[i].length){
				int idx = (int) u.random();
				int index = Arrays.binarySearch(rec_hotspot[i], idx);
				if (index < 0) {
					rec_hotspot[i][c++] = idx;
				}
			};
			Arrays.sort(rec_hotspot[i]);
			Tools.Fill_Matrix(rec_frac, i, rec_hotspot[i], 1);
			rec_frac[i][0] = 0.5;
		}
		
		if(DEBUG) print();
	}

	private void generateFounderGenotype(int[][] pair) {
		for (int i = 0; i < 2; i++) {
			for (int k = 0; k < N_snp; k++) {
				if (rnd.nextFloat() < snp_panel[k]) {
					pair[k][i] = 0;
				} else {
					pair[k][i] = 1;
				}
			}
		}
	}

	private void print() {
		if (DEBUG) {
			for (int i = 0; i < rec_hotspot.length; i++) {
				for (int j = 0; j < rec_hotspot[i].length; j++) {
					System.out.print(rec_hotspot[i][j] + " ");
				}
				System.out.println();
			}

			for (int i = 0; i < rec_frac.length; i++) {
				for (int j = 0; j < rec_frac[i].length; j++) {
					System.out.print(rec_frac[i][j] + " ");
				}
				System.out.println();
			}
		}
	}

	public static void main(String[] args) {
		double[] snp_freq = { 0.5, 0.5, 0.1, 0.1, 0.1 };
		HotSpot hs = new HotSpot(snp_freq);
		hs.GenerateChromosome(false);
	}
}
