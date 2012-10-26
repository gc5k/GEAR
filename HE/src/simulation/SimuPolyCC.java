package simulation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Calendar;

import parameter.Parameter;
import util.FileProcessor;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistributionImpl;
import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.random.RandomDataImpl;
import org.apache.commons.math.stat.StatUtils;

public class SimuPolyCC {

	private byte byte1 = 108;
	private byte byte2 = 27;
	private byte byte3 = 1;
	
	private NormalDistributionImpl norm;
	private RandomDataImpl rnd;
	private long seed = 2011;

	private int M = 1;
	private int M_null = 0;
	private boolean U = false;
	private int sample = 1000;
	private int N_case = 500;
	private int N_control = 500;

	private double[][] genotype;
	private double[] BV;
	private double[] phenotype;

	private double[] risk;

	private double[] freq;
	private double[] DPrime;
	private double ld;
	private double[] LD;

	private double vy = 1;
	private double h2 = 0.5;
	private double E = 0.5;
	private double K = 0.05;
	private String out = "Poly";

	private String A1 = "A";
	private String A2 = "C";
	public static StringBuilder LOG = new StringBuilder();

	public SimuPolyCC(Parameter P) {

		M = P.polyLoci;
		M_null = P.polyLociNull;
		U = P.polyU;
		ld = P.polyLD;
		seed = P.simuSeed;
		sample = P.simuCC[0] + P.simuCC[1];
		N_case = P.simuCC[0];
		N_control = P.simuCC[1];
		h2 = P.simuHsq;
		K = P.simuK;
		out = P.out;

		E = Math.sqrt(1 - h2);

		rnd = new RandomDataImpl();
		rnd.reSeed(seed);
		norm = new NormalDistributionImpl(0, vy);

		freq = new double[M];
		DPrime = new double[M - 1];

		Arrays.fill(freq, P.polyFreq);
		Arrays.fill(DPrime, ld);
		LD = CalculateDprime(freq, DPrime);

		genotype = new double[N_case + N_control][M];
		phenotype = new double[N_case + N_control];
		BV = new double[N_case + N_control];
		risk = new double[N_case + N_control];

		Calendar calendar = Calendar.getInstance();
		LOG.append("\nThe analysis was implemented at: " + calendar.getTime()
				+ "\n");
		LOG.append("seed: " + seed + "\n");
		LOG.append("MAF: " + P.polyFreq + "\n");
		LOG.append("Marker: " + M + "\n");
		LOG.append("Null Marker: " + M_null + "\n");
		if (Parameter.polyEffectFlag) {
			LOG.append("genetic effect file: " + P.polyEffectFile + "\n");
		} else {
			LOG.append("Uniform Effect: " + U + "\n");
		}
		LOG.append("Sample size: " + sample + "\n");
		LOG.append("case: " + N_case + "\n");
		LOG.append("Control: " + N_control + "\n");
		LOG.append("K: " + K + "\n");
		LOG.append("h2: " + h2 + "\n");
		LOG.append("out: " + out + "\n");
		if (Parameter.makebedFlag) {
			LOG.append("make bed");
		}
		LOG.append("\n");
	}

	public static void main(String[] args) {

		Parameter p = new Parameter();
		p.commandListenor(args);

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
			e.printStackTrace();
		}
		log.println(LOG.toString());
		log.close();

		System.out.println(LOG.toString());

	}

	public void GenerateSample() {
		GenerateSampleSelection();
		if(Parameter.makebedFlag) {
			writeBFile();
		} else {
			writeFile();
		}
		writeLog();
	}


	public void GenerateSampleSelection() {
		int count_case = 0;
		int count_control = 0;
		int count = 0;
		
		RealMatrix effect = null;
		if (Parameter.polyEffectFlag) {
			effect = readEffects();
		} else {
			effect = GenerateEffects();
		}

		norm = new NormalDistributionImpl(0, 1);

		int c = 0;
		while (count < sample) {
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

			if (1 - liability < K) {
				if (count_case < N_case) {
					BV[count_case] = bv;
					phenotype[count_case] = L;
					genotype[count_case] = chr.getColumn(0);
					risk[count_case] = liability;
					count_case++;
				} else {
					continue;
				}
			} else {
				if (count_control < N_control) {
					BV[N_case + count_control] = bv;
					phenotype[N_case + count_control] = L;
					genotype[N_case + count_control] = chr.getColumn(0);
					risk[N_case + count_control] = liability;
					count_control++;
				} else {
					continue;
				}
			}

			count++;
		}
		LOG.append("total individuals visited: " + c + "\n");
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
			for (int i = 0; i < effect.length - M_null; i++) {
				double sigma_b = Math.sqrt((vy * h2)
						/ ((M-M_null) * 2 * freq[i] * (1 - freq[i])));
				effect[i] = sigma_b;
			}
		} else {
			for (int i = 0; i < effect.length - M_null; i++) {
				double sigma_b = Math.sqrt((vy * h2)
						/ ((freq.length-M_null) * 2 * freq[i] * (1 - freq[i])));
				effect[i] = rnd.nextGaussian(0, sigma_b);
			}
		}
		
		if	(Parameter.simuOrderFlag) {
			Arrays.sort(effect);
		}

		PrintWriter eff = null;
		try {
			eff = new PrintWriter(new BufferedWriter(new FileWriter(out
					+ ".rnd")));
		} catch (IOException e) {
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
					v[i][j] = d < (f1 * f2 + LD[i - 1]) / f1 ? v[i - 1][j]
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

	public void writeBFile() {
		DataOutputStream bedout = null;
		PrintWriter fam = null;
		PrintWriter bim = null;
		PrintWriter phe = null;
		PrintWriter cov = null;
		PrintWriter geno = null;
		try {
			bedout = new DataOutputStream (new FileOutputStream(out
					+ ".bed"));
			fam = new PrintWriter(new BufferedWriter(new FileWriter(out + ".fam")));
			bim = new PrintWriter(new BufferedWriter(new FileWriter(out
					+ ".bim")));

			phe = new PrintWriter(new BufferedWriter(new FileWriter(out
					+ ".phe")));
			cov = new PrintWriter(new BufferedWriter(new FileWriter(out
					+ ".cov")));
			geno = new PrintWriter(new BufferedWriter(new FileWriter(out + ".add")));
		} catch (IOException e) {
			e.printStackTrace();
		}

		for (int i = 0; i < genotype.length; i++) {
			if (i < N_case) {
				fam.print("case_" + i + " ");
				fam.print(1 + " ");
				fam.print(0 + " ");
				fam.print(0 + " ");
				fam.print(1 + " ");

				fam.println(2 + " ");
			} else {
				fam.print("control_" + i + " ");
				fam.print(1 + " ");
				fam.print(0 + " ");
				fam.print(0 + " ");
				fam.print(1 + " ");

				fam.println(1 + " ");
			}

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

		for (int i = 0; i < genotype.length; i++) {
			for (int j = 0; j < genotype[i].length; j++) {
				geno.print(((int) genotype[i][j] + 1 ) + " "); 
			}
			geno.println();
		}

		geno.close();

		try {
			bedout.writeByte(byte1);
			bedout.writeByte(byte2);
			bedout.writeByte(byte3);
			for (int i = 0; i < M; i++) {
				byte gbyte = 0;
				int idx = 0;
				for (int j = 0; j < sample; j++) {
					int g = (int) genotype[j][i] + 1;
					switch(g) {
						case 0: g = 0; break;
						case 1: g = 2; break;
						case 2: g = 3; break;
						default: g = 1; break; //missing
					}

					g <<= 2 * idx;
					gbyte |= g;
					idx++;

					if (j != (sample - 1) ) {
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

		phe.close();
		bim.close();
		fam.close();
		cov.close();

	}

	public void writeFile() {
		PrintWriter pedout = null;
		PrintWriter map = null;
		PrintWriter phe = null;
		PrintWriter cov = null;
		PrintWriter geno = null;

		try {
			pedout = new PrintWriter(new BufferedWriter(new FileWriter(out
					+ ".ped")));
			map = new PrintWriter(new BufferedWriter(new FileWriter(out
					+ ".map")));
			phe = new PrintWriter(new BufferedWriter(new FileWriter(out
					+ ".phe")));
			cov = new PrintWriter(new BufferedWriter(new FileWriter(out
					+ ".cov")));
			geno = new PrintWriter(new BufferedWriter(new FileWriter(out + ".add")));
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

		for (int i = 0; i < genotype.length; i++) {
			for (int j = 0; j < genotype[i].length; j++) {
				geno.print(((int) genotype[i][j] + 1) + " "); 
			}
			geno.println();
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

	public double[] CalculateDprime(double[] f, double[] dprime) {
		double[] D = new double[dprime.length];

		for (int i = 0; i < D.length; i++) {
			if(dprime[i]>0) {
				D[i]=dprime[i]*Math.min(f[i]*(1-f[i+1]), f[i+1]*(1-f[i]));
			} else {
				D[i]=dprime[i]*Math.min(f[i]*f[i+1], (1-f[i])*(1-f[i+1]));
			}
		}

		return D;

	}
	
	public RealMatrix readEffects() {
		BufferedReader reader = FileProcessor.FileOpen(Parameter.polyEffectFile);
		double[] effect = new double[M];
		int c= 0;
		String line = null;
		try {
			while ((line = reader.readLine()) != null) {
				line.trim();
				String[] l = line.split(Parameter.whitespace);
				if(l.length < 1) continue;
				if( c < (M - M_null) ) {
					effect[c++] = Double.parseDouble(l[0]);
				}
			}
		} catch (IOException e) {
			e.printStackTrace(System.err);
			System.exit(0);
		}
		RealMatrix Eff = new Array2DRowRealMatrix(effect);
		
		if (Parameter.simuOrderFlag) {
			Arrays.sort(effect);
		}

		PrintWriter eff = null;
		try {
			eff = new PrintWriter(new BufferedWriter(new FileWriter(out
					+ ".rnd")));
		} catch (IOException e) {
			e.printStackTrace();
		}

		for (int i = 0; i < effect.length; i++) {
			eff.println("rs" + i + " " + A2 + " " + effect[i]);
		}
		eff.close();

		return Eff;
	}

}
