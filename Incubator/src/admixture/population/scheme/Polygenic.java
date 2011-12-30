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

public class Polygenic {

	private NormalDistributionImpl norm;
	private RandomDataImpl rnd;
	private long seed = 2011;

	private int M = 1;
	private int sample = 1000;
	private int N_case = 500;
	private int N_control = 500;

	private double[][] genotype;
	private double[] BV;
	private double[] phenotype;
	private double[] Liab;

	private double[] freq;
	private double[] LD;
	private double ld;
	private double[] Dprime;

	private double vy = 1;
	private double h2 = 0.5;
	private double E = 0.5;
	private double K = 0.05;

	private String out = "Poly";

	private String A1 = "A";
	private String A2 = "C";
	public static StringBuilder LOG = new StringBuilder();
	
	public Polygenic(PolygenicPar P) {

		M = P.marker;
		ld = P.ld;
		seed = P.seed;
		sample = P.sample;
		N_case = P.cs;
		h2 = P.h2;
		K = P.K;
		out = P.out;

		N_control = sample - N_case;
	    E = Math.sqrt(1-h2);

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
		Liab = new double[N_case + N_control];

		Calendar calendar = Calendar.getInstance();
		LOG.append("\nThe analysis was implemented at: " + calendar.getTime() + "\n");
		LOG.append("Marker: " + M + "\n");
		LOG.append("LD: " + ld + "\n");
		LOG.append("Sample size: " + sample + "\n");
		LOG.append("case: " + N_case + "\n");
		LOG.append("Control: " + N_control + "\n");
		LOG.append("K: " + K + "\n");
		LOG.append("h2: " + h2 + "\n");
		LOG.append("out: " + out + "\n");
		LOG.append("\n");
	}

	public static void main(String[] args) {

		PolygenicPar p = new PolygenicPar();
		p.commandListenor(args);

		Polygenic Poly = new Polygenic(p);
		Poly.GenerateSample();
		Poly.writeFile();
		Poly.writeLog();
		
	}
	
	public void writeLog() {

		Calendar calendar = Calendar.getInstance();
		LOG.append("\nThe analysis was completed at: " + calendar.getTime() + "\n");
		PrintWriter log = null;
		try {
			log = new PrintWriter(new BufferedWriter(new FileWriter(out + ".plog")));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		log.println(LOG.toString());
		log.close();
		
		System.out.println(LOG.toString());
		
	}

	public void GenerateSample() {

		int count_case = 0;
		int count_control = 0;
		int count = 0;
		RealMatrix effect = GenerateEffects();
		
		norm = new NormalDistributionImpl(0, 1);

		int c = 0;
		while(count < sample) {
			c++;
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

			if ( 1 - liability < K) {
				if (count_case < N_case) {
					BV[count_case] = bv;
					phenotype[count_case] = L;
					genotype[count_case] = chr.getColumn(0);
					Liab[count_case] = liability;
					count_case++;
				} else {
					continue;
				}
			} else {
				if (count_control < N_control) {
					BV[N_case + count_control] = bv;
					phenotype[N_case + count_control] = L;
					genotype[N_case + count_control] = chr.getColumn(0);
					Liab[N_case + count_control] = liability;
					count_control++;
				} else {
					continue;
				}
			}

			count++;
		}
		LOG.append("total individuals visited: " + c + "\n");

	}

	public RealMatrix GenerateEffects() {

		double[] effect = new double[freq.length];
		for (int i = 0; i < effect.length; i++) {
			double sigma_b = Math.sqrt((vy * h2)/ (freq.length * 2 * freq[i] * (1-freq[i])));
			effect[i] = rnd.nextGaussian(0, sigma_b);
		}

		PrintWriter eff = null;
		try {
			eff = new PrintWriter(new BufferedWriter(new FileWriter(out + ".rnd")));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		for(int i = 0; i < effect.length; i++) {
			eff.println("rs" + i + " " + A2 + " " + effect[i]);
		}
		eff.close();

		RealMatrix Eff = new Array2DRowRealMatrix(effect);
		
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
					int a = (int) v[i-1][j];
					double f1 = a == 0 ? freq[i-1] : (1-freq[i-1]);
					double f2 = a == 0 ? freq[i] : (1-freq[i]);
					v[i][j] = d < (f1 * f2 + Dprime[i-1]) / f1 ? v[i-1][j] : (1-v[i-1][j]);
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
			pedout = new PrintWriter(new BufferedWriter(new FileWriter(out + ".ped")));
			map = new PrintWriter(new BufferedWriter(new FileWriter(out + ".map")));
			phe = new PrintWriter(new BufferedWriter(new FileWriter(out + ".phe")));
			cov = new PrintWriter(new BufferedWriter(new FileWriter(out + ".cov")));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		for (int i = 0; i < genotype.length; i++) {
			if ( i < N_case) {
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

			for (int j = 0; j < genotype[i].length; j++ ) {
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

			if( i < N_case) {
				phe.print("case_" + i + " " + 1 + " " + 2 + " ");
				cov.print("case_" + i + " " + 1 + " " + 2 + " ");
			} else {
				phe.print("control_" + i + " " + 1 + " " + 1 + " ");
				cov.print("control_" + i + " " + 1 + " " + 1 + " ");
			}			
			phe.println();

			cov.print(Liab[i] + " ");
			cov.println(BV[i]);
		}

		for (int i = 0; i < M; i ++) {
			map.print(1 + " ");
			map.print("rs" + i + " ");
			map.print(i/(M * 1.0) + " ");
			map.println(i* 100);
		}
		
		pedout.close();
		phe.close();
		map.close();
		cov.close();
	}
	
	public double[] CalculateDprime(double[] f, double[] cor) {
		double[] D = new double[cor.length];

		for(int i = 0; i < D.length; i++) {
			D[i] = cor[i]*Math.sqrt(f[i] * (1-f[i]) * f[i+1] * (1-f[i+1]));
		}

		return D;

	}
}
