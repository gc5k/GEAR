package admixture.population.scheme;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Calendar;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistributionImpl;
import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.random.RandomDataImpl;
import org.apache.commons.math.stat.StatUtils;

public class PolygenicII {

	private NormalDistributionImpl norm;
	private RandomDataImpl rnd;
	private long seed = 2011;

	private int M = 1;
	private boolean U = false;
	private int sample = 1000;
	private int N_case = 500;
	private int N_control = 500;
	private int N_total;

	private double[][] genotype;
	private double[][] genotype_p;
	private double[] BV;
	private double[] BV_p;
	private double[] phenotype;
	private double[] phenotype_p;

	private double[] risk;

	private double[] freq;
	private double[] LD;
	private double ld;
	private double[] Dprime;

	private double vy = 1;
	private double h2 = 0.5;
	private double E = 0.5;
	private double K = 0.1;
	private boolean noselection = false;
	private String out = "Poly";

	private String A1 = "A";
	private String A2 = "C";
	public static StringBuilder LOG = new StringBuilder();

	public PolygenicII(NewPolygenicPar P) {

		M = P.marker;
		U = P.U;
		ld = P.ld;
		seed = P.seed;
		sample = P.sample;
		N_case = P.cs;
		h2 = P.h2;
		K = P.K;
		noselection = P.noselection;
		out = P.out;

		N_control = sample - N_case;
		E = Math.sqrt(1 - h2);

		rnd = new RandomDataImpl();
		rnd.reSeed(seed);
		norm = new NormalDistributionImpl(0, vy);

		freq = new double[M];
		LD = new double[M - 1];

		Arrays.fill(freq, 0.5);
		Arrays.fill(LD, ld);
		Dprime = CalculateDprime(freq, LD);

		genotype = new double[N_case + N_control][M];
		phenotype = new double[N_case + N_control];
		BV = new double[N_case + N_control];
		risk = new double[N_case + N_control];

		Calendar calendar = Calendar.getInstance();
		LOG.append("\nThe analysis was implemented at: " + calendar.getTime()
				+ "\n");
		LOG.append("seed: " + seed + "\n");
		LOG.append("Marker: " + M + "\n");
		LOG.append("Uniform Effet: " + U + "\n");
		LOG.append("LD: " + ld + "\n");
		LOG.append("Sample size: " + sample + "\n");
		LOG.append("case: " + N_case + "\n");
		LOG.append("Control: " + N_control + "\n");
		LOG.append("K: " + K + "\n");
		LOG.append("No selection: " + noselection + "\n");
		LOG.append("h2: " + h2 + "\n");
		LOG.append("out: " + out + "\n");
		LOG.append("\n");
	}

	public static void main(String[] args) {

		NewPolygenicPar p = new NewPolygenicPar();
		p.commandListenor(args);

		PolygenicII Poly = new PolygenicII(p);

		Poly.GenerateSample();
		Poly.writeFile();
		Poly.writeLog();

	}

	public void writeLog() {

		Calendar calendar = Calendar.getInstance();
		LOG.append("\nThe analysis was completed at: " + calendar.getTime()
				+ "\n");
		PrintWriter log = null;
		try {
			log = new PrintWriter(new BufferedWriter(new FileWriter(out
					+ ".plog")));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		log.println(LOG.toString());
		log.close();

		System.out.println(LOG.toString());

	}

	public void GenerateSample() {

		double total = Math.ceil(N_case / K);
		N_total = (new Double(total)).intValue();

		RealMatrix effect = GenerateEffects();
		BV_p = new double[N_total];
		genotype_p = new double[N_total][2];
		phenotype_p = new double[N_total];

		double BVsq = 0;
		double sum_bv = 0;
		for (int i = 0; i < N_total; i++ ) {
			RealMatrix chr = SampleChromosome();
			RealMatrix res = chr.transpose().multiply(effect);
			BV_p[i] = res.getEntry(0, 0);
			BVsq += BV_p[i] * BV_p[i];
			sum_bv += BV_p[i];
			genotype_p[i] = chr.getColumn(0);
		}
		double mean_bv = sum_bv/N_total;
		double sd_bv = Math.sqrt(BVsq/N_total - (mean_bv) * (mean_bv));

		for (int i = 0; i < N_total; i++) {
			BV_p[i] = Math.sqrt(h2) * (BV_p[i] - mean_bv) / sd_bv;
		}
		System.out.println("mean: " + BVsq + " " + mean_bv + " " + sd_bv);
		GenerateSampleSelection();
	}

	public void GenerateSampleSelection() {
		double[] pp = new double[N_total];
		for (int j = 0; j < N_total; j++) {
			phenotype_p[j] = BV_p[j] + rnd.nextGaussian(0, E);
			pp[j] = phenotype_p[j];
		}
		
		Arrays.sort(pp);
		double T = pp[N_total - 1 - N_case];
		
		int[] case_index = new int[N_case];
		int[] ctrl_index = new int[N_control];
		int idx_cs = 0;
		int idx_ctrl = 0;
		for (int j = 0; j < N_total; j++) {
			if(phenotype_p[j] >= T && idx_cs < N_case) {
				case_index[idx_cs++]= j;
			}
			if(phenotype_p[j] < T && idx_ctrl < N_control) {
				ctrl_index[idx_ctrl++] = j;
			}
		}
		
		for (int j = 0; j < N_case; j++) {
			System.out.println("cs " + j + ": " +BV_p[case_index[j]] + " " + phenotype_p[case_index[j]]);
		}
		for (int j = 0; j < N_control; j++) {
			System.out.println("ctrl " + j + ": " + BV_p[ctrl_index[j]] + " " + phenotype_p[ctrl_index[j]]);
		}
		
		System.out.println(N_total + " " + T);
	}

	public void GenerateSampleNoSelection() {

		int count_case = 0;
		int count_control = 0;
		int count = 0;
		RealMatrix effect = GenerateEffects();

		norm = new NormalDistributionImpl(0, 1);

		while (count < sample) {
			RealMatrix chr = SampleChromosome();
			RealMatrix res = chr.transpose().multiply(effect);

			double bv = res.getEntry(0, 0);
			double L = bv + rnd.nextGaussian(0, E);
			double liability = 0;

			try {
				liability = norm.cumulativeProbability(L);
			} catch (MathException e) {
				e.printStackTrace();
			}

			if (1 - liability < K) {
				BV[count_case] = bv;
				phenotype[count_case] = L;
				genotype[count_case] = chr.getColumn(0);
				risk[count_case] = liability;
				count_case++;
			} else {
				count_control++;
				BV[sample - count_control] = bv;
				phenotype[sample - count_control] = L;
				genotype[sample - count_control] = chr.getColumn(0);
				risk[sample - count_control] = liability;
			}

			count++;
		}
		N_case=count_case;
		LOG.append("total individuals visited (no selection): " + count + " (affected=" + N_case + ")\n");

	}

	public RealMatrix GenerateEffects() {

		double[] effect = new double[M];
		if (U) {
			for (int i = 0; i < effect.length; i++) {
				double sigma_b = Math.sqrt((vy * h2)
						/ (M * 2 * freq[i] * (1 - freq[i])));
				effect[i] = sigma_b;
			}
		} else {
			for (int i = 0; i < effect.length; i++) {
				double sigma_b = Math.sqrt((vy * h2)
						/ (freq.length * 2 * freq[i] * (1 - freq[i])));
				effect[i] = rnd.nextGaussian(0, sigma_b);
			}
		}

		PrintWriter eff = null;
		try {
			eff = new PrintWriter(new BufferedWriter(new FileWriter(out
					+ ".rnd")));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		for (int i = 0; i < effect.length; i++) {
			eff.println("rs" + i + " " + A2 + " " + effect[i]);
		}
		eff.close();

		RealMatrix Eff = new Array2DRowRealMatrix(effect);
		getMean(Eff);
		return Eff;

	}

	public RealMatrix SampleChromosome() {

		double[] g = new double[freq.length];
		double[][] v = new double[freq.length][2];
		for (int i = 0; i < freq.length; i++) {
			for (int j = 0; j < 2; j++) {
				double r = rnd.nextUniform(0, 1);
				double z = r < freq[i] ? 0 : 1;
				if (i == 0) {
					v[i][j] = z;
				} else {
					double d = rnd.nextUniform(0, 1);
					int a = (int) v[i - 1][j];
					double f1 = a == 0 ? freq[i - 1] : (1 - freq[i - 1]);
					double f2 = a == 0 ? freq[i] : (1 - freq[i]);
					v[i][j] = d < (f1 * f2 + Dprime[i - 1]) / f1 ? v[i - 1][j]
							: (1 - v[i - 1][j]);
				}
			}
			g[i] = v[i][0] + v[i][1] - 1;
		}

		RealMatrix chr = new Array2DRowRealMatrix(g);

		return chr;

	}

	public double getMean(RealMatrix effect) {

		int local_sample = 10000;
		double[] P = new double[local_sample];
		int c = 0;
		while(c < local_sample) {
			
			RealMatrix chr = SampleChromosome();
			RealMatrix res = chr.transpose().multiply(effect);
			P[c++] = res.getEntry(0, 0) + rnd.nextGaussian(0, E);

		}

		System.out.println(StatUtils.mean(P) + " " + StatUtils.variance(P));
		
		return StatUtils.mean(P);

	}

	public void writeFile() {
		PrintWriter pedout = null;
		PrintWriter map = null;
		PrintWriter phe = null;
		PrintWriter cov = null;
		try {
			pedout = new PrintWriter(new BufferedWriter(new FileWriter(out
					+ ".ped")));
			map = new PrintWriter(new BufferedWriter(new FileWriter(out
					+ ".map")));
			phe = new PrintWriter(new BufferedWriter(new FileWriter(out
					+ ".phe")));
			cov = new PrintWriter(new BufferedWriter(new FileWriter(out
					+ ".cov")));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		for (int i = 0; i < genotype.length; i++) {
			if (i < N_case) {
				pedout.print("case_" + i + " ");
				pedout.print(1 + " ");
				pedout.print(0 + " ");
				pedout.print(0 + " ");
				pedout.print(1 + " ");

				pedout.print(2 + " ");
			} else {
				pedout.print("control_" + i + " ");
				pedout.print(1 + " ");
				pedout.print(0 + " ");
				pedout.print(0 + " ");
				pedout.print(1 + " ");

				pedout.print(1 + " ");
			}

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

			if (i < N_case) {
				phe.print("case_" + i + " " + 1 + " " + 2 + " ");
				cov.print("case_" + i + " " + 1 + " " + 2 + " ");
			} else {
				phe.print("control_" + i + " " + 1 + " " + 1 + " ");
				cov.print("control_" + i + " " + 1 + " " + 1 + " ");
			}
			phe.println();

			cov.print(risk[i] + " ");
			cov.print(BV[i] + " ");
			cov.println(phenotype[i]);
		}

		for (int i = 0; i < M; i++) {
			map.print(1 + " ");
			map.print("rs" + i + " ");
			map.print(i / (M * 1.0) + " ");
			map.println(i * 100);
		}

		pedout.close();
		phe.close();
		map.close();
		cov.close();
	}

	public double[] CalculateDprime(double[] f, double[] cor) {
		double[] D = new double[cor.length];

		for (int i = 0; i < D.length; i++) {
			D[i] = cor[i]
					* Math.sqrt(f[i] * (1 - f[i]) * f[i + 1] * (1 - f[i + 1]));
		}

		return D;

	}
}
