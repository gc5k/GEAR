package gear.subcommands.simulationdipop;

import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Arrays;

import org.apache.commons.math.MathException;
import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.random.RandomDataImpl;
import org.apache.commons.math.stat.StatUtils;

import gear.ConstValues;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.Logger;
import gear.util.Sample;
import gear.util.pop.PopStat;

public class SimulationDiPopCommandImpl extends CommandImpl {

	@Override
	public void execute(CommandArguments cmdArgs) {
		dpArgs = (SimulationDiPopCommandArguments) cmdArgs;

		sample = dpArgs.getSampleSize();
		tsample = sample[0] + sample[1];

		h2 = dpArgs.getHsq();
		M = dpArgs.getTotalMarkerNum();
		Mset = dpArgs.getMarkerNum();
		
		fst = dpArgs.getFst();

		seed = dpArgs.getSeed();
		rnd.reSeed(seed);
		rep = dpArgs.getRep();

		getFreq();
		getEffect();
		getDPrime();
		calLD();
		generateSampleNoSelection();

		if (dpArgs.isMakeBed()) {
			writeBFile();
		} else {
			writeFile();
		}
		writeEffFile();
	}

	private void getFreq() {

		dFreq = new double[2][M];
		FREQ = new double[M];

		int[] mIdx = new int[fst.length];
		mIdx[0] = Mset[0];
		for (int i = 1; i < Mset.length; i++) {
			mIdx[i] = mIdx[i-1] + Mset[i];
		}

		for (int i = 0; i < M; i++) {
			double f = rnd.nextUniform(dpArgs.getFreqRangeLow(), dpArgs.getFreqRangeHigh());
			int k = 0;
			while (i > mIdx[k]) {
				k++;
			}
			try {
				dFreq[0][i] = rnd.nextBeta(f*(1-fst[k])/fst[k], (1-f)*(1-fst[k])/fst[k]);
				dFreq[1][i] = rnd.nextBeta(f*(1-fst[k])/fst[k], (1-f)*(1-fst[k])/fst[k]);
				FREQ[i] = sample[0]/tsample * dFreq[0][i] + sample[1]/tsample * dFreq[0][i];
			} catch (MathException e) {
				e.printStackTrace();
			}
		}
	}

	private void generateSampleNoSelection() {
		DecimalFormat fmt = new DecimalFormat("#.###E0");

		RealMatrix Meffect = new Array2DRowRealMatrix(effect);
		genotype = new double[tsample][M];
		phenotype = new double[tsample][rep];
		BV = new double[tsample];

		for (int i = 0; i < tsample; i++) {
			int idxFst = i < sample[0] ? 0:1;
			RealMatrix chr = SampleChromosome(idxFst);
			RealMatrix genoEff = chr.transpose().multiply(Meffect);

			double bv = genoEff.getEntry(0, 0);

			BV[i] = bv;
			genotype[i] = chr.getColumn(0);
		}

		double H2 = dpArgs.getHsq();
		if (H2 == 0) {
			Arrays.fill(BV, 0);
		}

		double vg = StatUtils.variance(BV);
		// rescale the phenotype to get the heritability and residual
		double ve = H2 == 0 ? 1 : vg * (1 - H2) / H2;
		double E = Math.sqrt(ve);
		Logger.printUserLog("Vg=" + fmt.format(vg));
		for (int i = 0; i < rep; i++) {
			double[] pv = new double[tsample];
			for (int j = 0; j < tsample; j++) {
				phenotype[j][i] = BV[j] + rnd.nextGaussian(0, E);
				pv[j] = phenotype[j][i];
			}
			double Vp = StatUtils.variance(pv);
			Logger.printUserLog("Vp=" + fmt.format(Vp) + "; hsq=" + fmt.format(vg / Vp) + " for replicate " + (i + 1));
		}
		Logger.printUserLog("Total individuals visited (no selection): " + BV.length + "\n");
	}

	public RealMatrix SampleChromosome(int fstIdx) {
		double[] g = new double[M];
		double[][] v = new double[M][2];
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < 2; j++) {
				double r = rnd.nextUniform(0, 1);
				double z = r < dFreq[fstIdx][i] ? 0 : 1;
				if (i == 0) {
					v[i][j] = z;
				} else {
					double d = rnd.nextUniform(0, 1);
					int a = (int) v[i - 1][j];
					double f1 = a == 0 ? dFreq[fstIdx][i - 1] : (1 - dFreq[fstIdx][i - 1]);
					double f2 = a == 0 ? dFreq[fstIdx][i] : (1 - dFreq[fstIdx][i]);
					v[i][j] = d < (f1 * f2 + LD[i - 1]) / f1 ? v[i - 1][j] : (1 - v[i - 1][j]);
				}
			}
			g[i] = v[i][0] + v[i][1] - 1;
		}

		RealMatrix chr = new Array2DRowRealMatrix(g);
		return chr;
	}

	private void getEffect() {
		effect = new double[M];

		Sample.setSeed(seed);

		for (int i = 0; i < M; i++) {
			effect[i] = rnd.nextGaussian(0, 1);
		}

		if (h2 == 0) {
			Arrays.fill(effect, 0);
		}
	}

	private void getDPrime() {
		dprime = new double[M - 1];
		if (dpArgs.isPlainLD()) {
			Arrays.fill(dprime, dpArgs.getLD());
		} else if (dpArgs.isRandLD()) {
			for (int i = 0; i < dprime.length; i++) {
				dprime[i] = rnd.nextUniform(-1, 1);
			}
		}
	}

	public void calLD() {
		LD = PopStat.CalcLDfromDPrime(FREQ, dprime);
	}

	public void writeBFile() {
		DataOutputStream bedout = null;
		PrintWriter fam = null;
		PrintWriter bim = null;
		PrintWriter phe = null;
		PrintWriter geno = null;
		PrintWriter breed = null;

		try {
			bedout = new DataOutputStream(new FileOutputStream(dpArgs.getOutRoot() + ".bed"));

			fam = new PrintWriter(new BufferedWriter(new FileWriter(dpArgs.getOutRoot() + ".fam")));
			bim = new PrintWriter(new BufferedWriter(new FileWriter(dpArgs.getOutRoot() + ".bim")));
			phe = new PrintWriter(new BufferedWriter(new FileWriter(dpArgs.getOutRoot() + ".phe")));
			geno = new PrintWriter(new BufferedWriter(new FileWriter(dpArgs.getOutRoot() + ".add")));
			breed = new PrintWriter(new BufferedWriter(new FileWriter(dpArgs.getOutRoot() + ".breed")));
		} catch (IOException e) {
			e.printStackTrace();
		}

		for (int i = 0; i < genotype.length; i++) {
			fam.print("sample_" + i + " ");

			fam.print(1 + " ");
			fam.print(0 + " ");
			fam.print(0 + " ");
			fam.print(1 + " ");

			fam.println(phenotype[i][0] + " ");

			phe.print("sample_" + i + " " + 1);
			breed.println("sample_" + i + " " + 1 + " " + BV[i]);

			for (int j = 0; j < rep; j++) {
				phe.print(" " + phenotype[i][j]);
			}
			phe.println();
		}

		for (int i = 0; i < genotype.length; i++) {
			for (int j = 0; j < genotype[i].length; j++) {
				geno.print(((int) genotype[i][j] + 1) + " ");
			}
			geno.println();
		}

		try {
			bedout.writeByte(ConstValues.PLINK_BED_BYTE1);
			bedout.writeByte(ConstValues.PLINK_BED_BYTE2);
			bedout.writeByte(ConstValues.PLINK_BED_BYTE3);
			for (int i = 0; i < M; i++) {
				byte gbyte = 0;
				int idx = 0;
				for (int j = 0; j < tsample; j++) {
					int g = (int) genotype[j][i] + 1;
					switch (g) {
					case 0:
						g = 0;
						break;
					case 1:
						g = 2;
						break;
					case 2:
						g = 3;
						break;
					default:
						g = 1;
						break; // missing
					}

					g <<= 2 * idx;
					gbyte |= g;
					idx++;

					if (j != (tsample - 1)) {
						if (idx == 4) {
							bedout.writeByte(gbyte);
							gbyte = 0;
							idx = 0;
						}
					} else {
						bedout.writeByte(gbyte);
					}
				}
			}
			bedout.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		for (int i = 0; i < M; i++) {
			bim.print(1 + " ");
			bim.print("rs" + i + " ");
			bim.print(i / (M * 1.0) + " ");
			bim.print(i * 100 + " ");
			bim.println(A1 + " " + A2);
		}

		geno.close();
		phe.close();
		bim.close();
		fam.close();
		breed.close();
	}

	public void writeFile() {
		PrintWriter pedout = null;
		PrintWriter map = null;
		PrintWriter phe = null;
		PrintWriter geno = null;
		PrintWriter breed = null;
		try {
			pedout = new PrintWriter(new BufferedWriter(new FileWriter(dpArgs.getOutRoot() + ".ped")));
			map = new PrintWriter(new BufferedWriter(new FileWriter(dpArgs.getOutRoot() + ".map")));
			phe = new PrintWriter(new BufferedWriter(new FileWriter(dpArgs.getOutRoot() + ".phe")));
			geno = new PrintWriter(new BufferedWriter(new FileWriter(dpArgs.getOutRoot() + ".add")));
			breed = new PrintWriter(new BufferedWriter(new FileWriter(dpArgs.getOutRoot() + ".breed")));
		} catch (IOException e) {
			Logger.handleException(e, "An exception occurred when writing files.");
		}

		for (int i = 0; i < genotype.length; i++) {
			pedout.print("sample_" + i + " ");
			pedout.print(1 + " ");
			pedout.print(0 + " ");
			pedout.print(0 + " ");
			pedout.print(1 + " ");
			pedout.print(phenotype[i][0] + " ");

			for (int j = 0; j < genotype[i].length; j++) {
				int g = (int) genotype[i][j];
				if (g == -1) {
					pedout.print(A1 + " " + A1 + "  ");
				} else if (g == 0) {
					pedout.print(A1 + " " + A2 + "  ");
				} else {
					pedout.print(A2 + " " + A2 + "  ");
				}
			}
			pedout.println();

			phe.print("sample_" + i + " " + 1);
			breed.println("sample_" + i + " " + 1 + " " + BV[i]);
			for (int j = 0; j < rep; j++) {
				phe.print(" " + phenotype[i][j]);
			}
			phe.println();
		}

		for (int i = 0; i < M; i++) {
			map.print(1 + " ");
			map.print("rs" + i + " ");
			map.print(i / (M * 1.0) + " ");
			map.println(i * 100);
		}

		for (int i = 0; i < genotype.length; i++) {
			for (int j = 0; j < genotype[i].length; j++) {
				geno.print(((int) genotype[i][j] + 1) + " ");
			}
			geno.println();
		}

		pedout.close();
		map.close();
		phe.close();
		geno.close();
		breed.close();
	}

	private void writeEffFile() {
		PrintWriter eff = null;
		try {
			eff = new PrintWriter(new BufferedWriter(new FileWriter(dpArgs.getOutRoot() + ".rnd")));
		} catch (IOException e) {
			Logger.handleException(e, "An exception occurred when writing files.");
		}

		for (int i = 0; i < M; i++) {
			eff.println("rs" + i + " " + A1 + " " + effect[i]);
		}
		eff.close();
	}

	private SimulationDiPopCommandArguments dpArgs;

	private RandomDataImpl rnd = new RandomDataImpl();
	private long seed;

	private int M;
	private int[] Mset;
	private double[] fst;
	private double[] FREQ;
	private double[][] dFreq;

	private int tsample;
	private int[] sample;
	private int rep;

	private double[][] genotype;
	private double[] BV;
	private double[][] phenotype;

	private double[] effect;

	private double[] dprime;
	private double[] LD;

	private double h2;

	private String A1 = "A";
	private String A2 = "C";
}
