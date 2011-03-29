package admixture.population.genome.chromosome;

import java.util.Random;

import admixture.AdmixtureConstant;
import admixture.population.genome.HotSpot;

/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */

public class ChromosomeGenerator {
	protected boolean DEBUG = false;
	protected int N_snp;
	protected Random rnd;
	protected double[] snp_panel;
	protected double[][] ancestry_snp_panel;
	protected double[] ancestry;
	protected double[] LD; // reserve for the future
	protected HotSpot hs;
	protected int[][] hotspot; // there two kinds of recombination fractions

	// 1: free of recombination that each element equals 0.5, hotspot=[0.5,0.5,0.5,...]
	// it is generated in the method recombinationFree()
	// 2: cross over between snps: hotspot looks like=[0, 0, 0, 1, 0, 0, 1, ...]
	// it is generated in the method recombination()

	// Note: hotspot[0] for paternal recombination, hotspot[1] for maternal recombination;

	// protected int[][][] parent_ancestry;
	// protected int[][][] offspring_ancestry;

	public ChromosomeGenerator(double[][] sp, double[] anc) {
		ancestry_snp_panel = sp;
		ancestry = new double[anc.length];
		double[] a = new double[anc.length];
		double s = 0;
		for (int i = 0; i < anc.length; i++) {
			s += anc[i];
			if (i == 0) {
				a[i] = anc[i] + 0;
			} else {
				a[i] += anc[i] + a[i - 1];
			}
		}

		for (int i = 0; i < anc.length - 1; i++) {
			ancestry[i] = a[i] / s;
		}
		ancestry[ancestry.length - 1] = 1;

		rnd = new Random(2011);
		hotspot = new int[2][];
		N_snp = ancestry_snp_panel[0].length;
	}

	public FamilySingleChromosome generateFamilySingleChromosome(int cID, int k, HotSpot h, boolean disease_linked) {
		hs = h;
		// parent_ancestry = new int[2][2][N_snp];
		int[][][] pg = generateFounderChr(AdmixtureConstant.Without_LD);

		// offspring_ancestry = new int[k][2][N_snp];
		int[][][] og = generateOffspringChr(pg, k);

		FamilySingleChromosome fsc = new FamilySingleChromosome(cID, pg, og, disease_linked, ancestry.length);
		return fsc;
	}

	public FamilySingleChromosome generateFamilySingleChromosome(int cID, int[][] p_g, int[][] m_g, int k, HotSpot h,
			boolean disease_linked) {
		hs = h;
		// parent_ancestry = new int[2][][];
		// parent_ancestry[0] = p_a;
		// parent_ancestry[1] = m_a;
		int[][][] pg = new int[2][][];
		pg[0] = p_g;
		pg[1] = m_g;

		// offspring_ancestry = new int[k][2][N_snp];
		int[][][] og = generateOffspringChr(pg, k);

		FamilySingleChromosome fsc = new FamilySingleChromosome(cID, pg, og, disease_linked, ancestry.length);
		return fsc;
	}

	public FamilySingleChromosome generateFamilySingleChromosomeII(int cID, int k, int[] father_hotspot,
			int[] mother_hotspot, int[][][] parent_anc, int[][][] p_g, boolean disease_linked) {
		// parent_ancestry = parent_anc;
		int[][][] pg = p_g;

		// offspring_ancestry = new int[k][2][N_snp];
		int[][][] og = generateOffspringChr(pg, k);

		FamilySingleChromosome fsc = new FamilySingleChromosome(cID, pg, og, disease_linked, ancestry.length);
		return fsc;
	}

	private int[][][] generateFounderChr(boolean ld) {
		int[][][] diploid = new int[2][2][N_snp];
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				double r = rnd.nextFloat();
				int idx = 0;
				while (r >= ancestry[idx] && idx < ancestry.length)
					idx++;
				for (int k = 0; k < N_snp; k++) {
					if (!ld) {
						diploid[i][j][k] = rnd.nextFloat() < ancestry_snp_panel[idx][k] ? idx * 2 + 0 : idx * 2 + 1;
					} else {
						// when there is LD pattern;
					}
				}
				// Arrays.fill(parent_ancestry[i][j], idx);
			}
		}
		return diploid;
	}

	private int[][][] generateOffspringChr(int[][][] pg, int N_kid) {// pg[0] for paternal chromosome, pg[1] for maternal.
		int[][][] diploid = new int[N_kid][2][N_snp];
		for (int o = 0; o < N_kid; o++) {
			hs.GenerateRecombination(AdmixtureConstant.free_recombination);
			hotspot[0] = hs.getHotSpot();
			hs.GenerateRecombination(AdmixtureConstant.free_recombination);
			hotspot[1] = hs.getHotSpot();
			for (int i = 0; i < hotspot.length; ++i) {
				int chromatid = rnd.nextBoolean() ? 0 : 1;
				for (int j = 0; j < hotspot[i].length - 1; j++) {
					chromatid = 1 - chromatid;
					for (int k = hotspot[i][j]; k <= hotspot[i][j + 1]; k++) {
						diploid[o][i][k] = pg[i][chromatid][k];
						// offspring_ancestry[o][i][k] = parent_ancestry[i][chromatid][k];
					}
				}
			}
		}
		return diploid;
	}

	public void setSeed(long s) {
		rnd.setSeed(s);
	}

	public static void main(String[] args) {
		double[][] snp_freq = { { 0.5, 0.5, 0.1, 0.1, 0.1 }, { 0.1, 0.1, 0.1, 0.5, 0.5 } };
		double[] anc = { 0.5, 0.5 };

		HotSpot hs = new HotSpot();
		hs.rev(snp_freq[0].length);
		hs.GenerateRecombination(AdmixtureConstant.free_recombination);
		int[] fh = hs.getHotSpot();
		hs.GenerateRecombination(AdmixtureConstant.free_recombination);
		int[] mh = hs.getHotSpot();

		ChromosomeGenerator CG = new ChromosomeGenerator(snp_freq, anc);
		CG.generateFamilySingleChromosome(0, 2, hs, true);
		hs.GenerateRecombination(AdmixtureConstant.free_recombination);
		fh = hs.getHotSpot();
		hs.GenerateRecombination(AdmixtureConstant.free_recombination);
		mh = hs.getHotSpot();
		CG.generateFamilySingleChromosome(0, 2, hs, true);

	}
}
