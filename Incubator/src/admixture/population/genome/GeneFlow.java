package admixture.population.genome;

import java.util.Arrays;
import java.util.Random;

import admixture.AdmixtureConstant;
import arsenal.Sample;

import jsc.combinatorics.Permutation;
import jsc.combinatorics.Permutations;

public class GeneFlow {
	private int[][][] pool;
	private int[][][] pool_ancestry;
	private int[][][] f_pool;
	private int[][][] f_pool_ancestry;

	private double[][] ancestry_snp_panel;
	int N_snp;
	private boolean ld = false;
	private Random rnd = new Random(2011);
	private double[] pop_prop;
	private int founder_size;

	private HotSpot hs;
	private int[][] hotspot = new int[2][];

	public GeneFlow(int ps, double[][] sp, double[] pp, HotSpot h) {
		founder_size = ps;
		N_snp = sp[0].length;
		ancestry_snp_panel = sp;
		pop_prop = new double[pp.length];
		System.arraycopy(pp, 0, pop_prop, 0, pp.length);
		for (int i = 0; i < pp.length - 1; i++) {
			pop_prop[i + 1] = pp[i + 1] + pop_prop[i];
		}
		hs = h;
		hs.rev(N_snp);
	}

	public int[][] generateAnFounder(int idx) {
		int[][] g = new int[2][N_snp];
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < N_snp; j++) {
				g[i][j] = rnd.nextFloat() < ancestry_snp_panel[idx][j] ? 0 : 1;
			}
		}
		return g;
	}

	public void founder() {
		pool = new int[founder_size][2][N_snp];
		pool_ancestry = new int[founder_size][2][N_snp];
		for (int i = 0; i < founder_size; i++) {
			for (int j = 0; j < 2; j++) {
				double r = rnd.nextFloat();
				int idx = 0;
				while (r >= pop_prop[idx] && idx < pop_prop.length)
					idx++;
				for (int k = 0; k < N_snp; k++) {
					if (!ld) {
						pool[i][j][k] = rnd.nextFloat() < ancestry_snp_panel[idx][k] ? 0 : 1;
					} else {
						// when there is LD pattern;
					}
				}
				Arrays.fill(pool_ancestry[i][j], idx);
			}
		}
	}

	private void MakePool() {

		int EA_size = EA_size();
		int AA_size = founder_size - EA_size;
		int[] AA_idx = Sample.SampleIndex(0, pool.length - 1, AA_size, false);

		f_pool = new int[pool.length][2][N_snp];
		f_pool_ancestry = new int[pool.length][2][N_snp];

		for (int i = 0; i < AA_idx.length; i++) {
			for (int j = 0; j < 2; j++) {
				System.arraycopy(pool[AA_idx[i]][j], 0, f_pool[i][j], 0, pool[AA_idx[i]][j].length);
				System.arraycopy(pool_ancestry[i][j], 0, f_pool_ancestry[i][j], 0, pool_ancestry[i][j].length);
			}
		}

		for (int i = 0; i < EA_size; i++) {
			int[][] g = generateAnFounder(1);
			System.arraycopy(g[0], 0, f_pool[i + AA_size][0], 0, g[0].length);
			System.arraycopy(g[1], 0, f_pool[i + AA_size][1], 0, g[1].length);
			Arrays.fill(f_pool[i + AA_size][0], 1);
			Arrays.fill(f_pool[i + AA_size][1], 1);
		}
	}

	public void mating() {
		MakePool();
		int[] idx = new int[f_pool.length];
		for (int i = 0; i < idx.length; i++) {
			idx[i] = i + 1;
		}
		Permutations ps = new Permutations(founder_size);
		Permutation p = ps.randomPermutation();
		int[] p_idx = p.toIntArray();
		for(int i = 0; i < p_idx.length; i++) p_idx[i]--;

		int[][][] pool = new int[f_pool.length][2][N_snp];
		int count = 0;
		for (int i = 0; i < f_pool.length / 2; i++) {
			int idxf = p_idx[i*2];
			int idxm = p_idx[i*2+1];
			for(int j = 0; j < 2; j++) {
				generateOffspringChr(count++, f_pool[idxf], f_pool[idxm], f_pool_ancestry[idxf], f_pool_ancestry[idxm]);
			}
		}
	}

	private void generateOffspringChr(int idx, int[][] fg, int[][] mg, int[][] fa, int[][] ma) {// pg[0] for paternal chromosome, pg[1] for
																									// maternal.
		int[][][] pg = new int[2][][];
		pg[0] = fg;
		pg[1] = mg;
		int[][][] pa = new int[2][][];
		pa[0] = fa;
		pa[1] = ma;

		hs.GenerateRecombination(AdmixtureConstant.free_recombination);
		hotspot[0] = hs.getHotSpot();
		hs.GenerateRecombination(AdmixtureConstant.free_recombination);
		hotspot[1] = hs.getHotSpot();
		for (int i = 0; i < hotspot.length; ++i) {
			int chromatid = rnd.nextBoolean() ? 0 : 1;
			for (int j = 0; j < hotspot[i].length - 1; j++) {
				chromatid = 1 - chromatid;
				for (int k = hotspot[i][j]; k <= hotspot[i][j + 1]; k++) {
					pool[idx][i][k] = pg[i][chromatid][k];
					pool_ancestry[idx][i][k] = pa[i][chromatid][k];
				}
			}
		}
	}

	public double[][] averageAncestry() {
		double[][] a = new double[founder_size][pop_prop.length];
		for(int i = 0; i < founder_size; i++) {
			for(int j = 0; j < N_snp; j++) {
				int org1 = pool_ancestry[i][0][j];
				int org2 = pool_ancestry[i][1][j];
				a[i][org1] += 0.5;
				a[i][org2] += 0.5;
			}
			for(int j = 0; j < a[i].length; j++) {
				a[i][j] /= N_snp;
			}
		}
		return a;
	}

	private int PopulationSizeNextGeneration() {
		return founder_size;
	}

	private int EA_size() {
		return (int) (founder_size * (pop_prop[1] - pop_prop[0]));
	}

	public static void main(String[] args) {
		double[][] snp = new double[2][100];

		Arrays.fill(snp[0], 0.1);
		Arrays.fill(snp[1], 0.8);

		double[] r = { 0.8, 0.2 };
		int ps = 100;
		HotSpot hs = new HotSpot();

		GeneFlow gf = new GeneFlow(ps, snp, r, hs);
		gf.founder();
//		gf.mating();
//		gf.mating();
		double[][] a = gf.averageAncestry();
		for(int i = 0; i < a.length; i++) {
			for(int j = 0; j < a[i].length; j++) {
				System.out.print(a[i][j] + " ");
			}
			System.out.println();
		}
	}
}
