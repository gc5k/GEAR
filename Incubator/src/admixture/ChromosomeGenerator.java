package admixture;

import java.util.ArrayList;
import java.util.Random;

public class ChromosomeGenerator {
	private boolean DEBUG = true;
	private int N_snp;
	private Random rnd = new Random(2011);
	private double[] snp_panel;
	private double[][][] post_prob;
	private double[] LD; // reserve for the future
	private double[][] rec_frac; // there two kinds of recombination fractions
	// 1: free of recombination that each element equals 0.5, rec_frac=[0.5,0.5,0.5,...]
	// it is generated in the method recombinationFree()
	// 2: cross over between snps: rec_frac looks like=[0, 0, 0, 1, 0, 0, 1, ...]
	// it is generated in the method recombination()
	
	// Note: rec_frac[0] for paternal recombination, rec_frac[1] for maternal recombination;

	private ArrayList<Family> habitat = new ArrayList<Family>();
	public static class Family {
		protected int[][][] p_g; // parents chromosomes;
								 // p_g[0] for dad, p_g[1] for mom
		protected int[][][] o_g; // offspring chromosomes;

		public Family(int[][][]p, int[][][] o) {
			p_g = p;
			o_g = o;
		}
	}

	public ChromosomeGenerator(double[] sp, double[][][] pp) {
		snp_panel = sp;
		post_prob = pp;

		N_snp = sp.length;
	}

	public void generatePedigree(int k, double[][] rf) {
		rec_frac = rf;
		int[][][] pg = new int[2][][];
		pg[0] = generateFounderChr(false);
		pg[1] = generateFounderChr(false);

		int[][][] og = new int[k][][];
		for (int i = 0; i < k; i++) {
			og[i] = generateOffspringChr(pg);
		}
		Family f = new Family(pg, og);
		habitat.add(f);
	}

	private int[][] generateFounderChr(boolean ld) {
		int[][] diploid = new int[2][N_snp];
		for (int i = 0; i < 2; i++) {
			for (int k = 0; k < N_snp; k++) {
				if(!ld) {
					diploid[k][i] = rnd.nextFloat() < snp_panel[k] ? 0:1;
				} else {
					// when there is LD pattern;
				}
			}
		}
		return diploid;
	}

	private int[][] generateOffspringChr(int[][][] pg) {//pg[0] for paternal chromosome, pg[1] for maternal.
		int[][] diploid = new int[2][N_snp];
		for (int i = 0; i < 2; ++i) {
			int chromatid = 0;
			for (int k = 0; k < N_snp; ++k) {
				double rd = rnd.nextFloat();
				if (rd < rec_frac[i][k]) {
					chromatid = 1 - chromatid;
				}
				diploid[i][k] = pg[i][chromatid][k];
			}
		}
		return diploid;
	}

	public static void main(String[] args) {
		double[] snp_freq = { 0.5, 0.5, 0.1, 0.1, 0.1 };
		HotSpot hs = new HotSpot(snp_freq);
		hs.GenerateChromosome(false);
	}
}
