package admixture.population.genome;

import java.util.Arrays;
import java.util.Random;

import admixture.AdmixtureConstant;
import admixture.population.AlleleFrequencyReader;
import arsenal.Sample;

public class GeneFlowByteFormat {
	private byte[][][] pool;
	private byte[][][] f_pool;

	private double[][] ancestry_snp_panel;
	int N_snp;
	int n_dim;
	private boolean ld = false;
	private Random rnd = new Random(2011);
	private double[] pop_prop;
	private int founder_size;

	private HotSpot hs;
	private int[][] hotspot = new int[2][];

	private int[] randomIndex;

	public GeneFlowByteFormat(int ps, double[][] sp, double[] pp, HotSpot h) {
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

		int n_dim = N_snp / 4;
		if (N_snp % 4 > 0) {
			n_dim++;
		}
		pool = new byte[founder_size][2][n_dim];
		f_pool = new byte[founder_size][2][n_dim];
	}

	public GeneFlowByteFormat(AlleleFrequencyReader afr, int ps, double[] pp, HotSpot h) {
		founder_size = ps;
		N_snp = afr.getNumberSNP();
		ancestry_snp_panel = afr.getAlleleFreq();

		pop_prop = new double[pp.length];
		System.arraycopy(pp, 0, pop_prop, 0, pp.length);
		for (int i = 0; i < pp.length - 1; i++) {
			pop_prop[i + 1] = pp[i + 1] + pop_prop[i];
		}

		hs = h;
		hs.rev(N_snp);

		int n = N_snp / 4;
		if (N_snp % 4 > 0) {
			n++;
		}
		pool = new byte[founder_size][2][n];
		f_pool = new byte[founder_size][2][n];
	}

	public void generateAnFounder(int idx, int slot) {
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < N_snp; j++) {
				byte g = (byte) (rnd.nextFloat() < ancestry_snp_panel[idx][j] ? (idx * 2 + 0) : (idx * 2 + 1));
				set_f_pool(slot, i, j, g);
			}
			// Arrays.fill(f_pool_ancestry[slot][0], 1);
		}
	}

	private void set_f_pool(int slot, int haploid, int j, byte g) {
		int loc = j / 4;
		int shift = (j % 4) * 2;
		f_pool[slot][haploid][loc] += g << shift;
	}

	private byte get_f_pool(int slot, int haploid, int j) {
		int loc = j / 4;
		int m_idx = j % 4;
		int shift = m_idx * 2;
		byte g = (byte) ((f_pool[slot][haploid][loc] >> shift) & 3);
		return g;
	}

	private void set_pool(int slot, int haploid, int j, byte g) {
		int loc = j / 4;
		int shift = (j % 4) * 2;
		f_pool[slot][haploid][loc] += g << shift;
	}

	private byte get_pool(int slot, int haploid, int j) {
		int loc = j / 4;
		int m_idx = j % 4;
		int shift = m_idx * 2;
		byte g = (byte) (pool[slot][haploid][loc] >> shift & 3);
		return g;
	}

	private void GenerateFounder() {

		for (int i = 0; i < founder_size; i++) {
			for (int j = 0; j < 2; j++) {
				double r = rnd.nextFloat();
				int idx = 0;
				while (r >= pop_prop[idx] && idx < pop_prop.length)
					idx++;
				for (int k = 0; k < N_snp; k++) {
					if (!ld) {
						byte g = (byte) (rnd.nextFloat() < ancestry_snp_panel[idx][k] ? (idx * 2 + 0) : (idx * 2 + 1));
						set_pool(i, j, k, g);
					} else {
						// when there is LD pattern;
					}
				}
				// Arrays.fill(pool_ancestry[i][j], idx);
			}
		}
	}

	private void MakePool() {
		int EA_size = EA_size();
		int AA_size = founder_size - EA_size;
		int[] AA_idx = Sample.SampleIndex(0, pool.length - 1, AA_size);

		for (int i = 0; i < AA_idx.length; i++) {
			for (int j = 0; j < 2; j++) {
				System.arraycopy(pool[AA_idx[i]][j], 0, f_pool[i][j], 0, pool[AA_idx[i]][j].length);
				// System.arraycopy(pool_ancestry[i][j], 0,
				// f_pool_ancestry[i][j], 0, pool_ancestry[i][j].length);
			}
		}

		for (int i = 0; i < EA_size; i++) {
			generateAnFounder(1, i + AA_size);
		}
	}

	public void mating(int Round) {
		GenerateFounder();
		for (int round = 0; round < Round - 1; round++) {
			MakePool();
			int[] p_idx = Sample.SampleIndex(0, founder_size - 1, founder_size);

			int count = 0;
			for (int i = 0; i < f_pool.length / 2; i++) {
				for (int j = 0; j < 2; j++) {
					generateOffspringChr(count++, p_idx[i * 2], p_idx[i * 2 + 1]);
				}
			}
		}
		MakePool();
	}

	private void generateOffspringChr(int idx, int idxf, int idxm) {
		// maternal.
		int[] IDX = { idxf, idxm };
		hs.GenerateRecombination(AdmixtureConstant.free_recombination);
		hotspot[0] = hs.getHotSpot();
		hs.GenerateRecombination(AdmixtureConstant.free_recombination);
		hotspot[1] = hs.getHotSpot();
		for (int i = 0; i < hotspot.length; ++i) {
			int chromatid = rnd.nextBoolean() ? 0 : 1;
			for (int j = 0; j < hotspot[i].length - 1; j++) {
				chromatid = 1 - chromatid;
				for (int k = hotspot[i][j]; k < hotspot[i][j + 1]; k++) {
					int loc = k / 4;
					pool[idx][i][loc] += get_f_pool(IDX[i], i, hotspot[i][j]);
				}
			}
		}
	}

	public double[][] averageAncestry() {
		double[][] a = new double[founder_size][pop_prop.length];
		for (int i = 0; i < founder_size; i++) {
			for (int j = 0; j < N_snp; j++) {
				byte g1 = get_pool(i, 0, j);
				byte g2 = get_pool(i, 1, j);
				int org1 = (int) (g1 / 2);
				int org2 = (int) (g2 / 2);
				a[i][org1] += 0.5;
				a[i][org2] += 0.5;
			}
			for (int j = 0; j < a[i].length; j++) {
				a[i][j] /= N_snp;
			}
		}
		return a;
	}

	public int NumberOfSNP() {
		return N_snp;
	}

	@SuppressWarnings("unused")
	private int PopulationSizeNextGeneration() {
		return founder_size;
	}

	private int EA_size() {
		return (int) (founder_size * (pop_prop[1] - pop_prop[0]));
	}

	public byte[][] getAnIndividualInPool(int idx) {
		if (idx < randomIndex.length) {
			return f_pool[randomIndex[idx]];
		} else {
			return null;
		}
	}

	public byte[][] sampleAnFounder(int pi) {
		byte[][] g = new byte[2][N_snp];
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < N_snp; j++) {
				int loc = j / 4;
				int shift = (j % 4) * 2;
				byte geno = (byte) (rnd.nextFloat() < ancestry_snp_panel[pi][j] ? (pi * 2 + 0) : (pi * 2 + 1));
				g[i][loc] += geno << shift;
			}
		}
		return g;
	}

	public void randomization() {
		randomIndex = Sample.SampleIndex(0, founder_size - 1, founder_size);
	}

	public static void main(String[] args) {
		double[][] snp = new double[2][1000];

		Arrays.fill(snp[0], 0.1);
		Arrays.fill(snp[1], 0.8);

		double[] r = { 0.975, 0.025 };
		int ps = 1000;
		HotSpot hs = new HotSpot();

		GeneFlowByteFormat gf = new GeneFlowByteFormat(ps, snp, r, hs);
		gf.mating(100);
		double[][] a = gf.averageAncestry();
		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < a[i].length; j++) {
				System.out.print(a[i][j] + " ");
			}
			System.out.println();
		}
	}
}
