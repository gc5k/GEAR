package simulation;

import gear.CmdArgs;
import gear.ConstValues;
import gear.util.FileUtil;
import gear.util.Logger;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Arrays;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.random.RandomDataImpl;
import org.apache.commons.math.stat.StatUtils;

public class SimuPolyQT
{

	private byte byte1 = 108;
	private byte byte2 = 27;
	private byte byte3 = 1;

	private RandomDataImpl rnd;
	private long seed = 2011;

	private int M = 100;
	private int M_null = 0;
	private boolean U = false;
	private int sample = 1000;

	private double[][] genotype;
	private double[] BV;
	private double[] phenotype;

	private double[] freq;
	private double[] DPrime;
	private double ld;
	private double[] LD;

	private double vy = 1;
	private double h2 = 0.5;
	private double E = 0.5;


	private String A1 = "A";
	private String A2 = "C";

	public SimuPolyQT()
	{
		M = CmdArgs.INSTANCE.polyLoci;
		M_null = CmdArgs.INSTANCE.polyLociNull;
		U = CmdArgs.INSTANCE.polyU;
		ld = CmdArgs.INSTANCE.polyLD;
		seed = CmdArgs.INSTANCE.simuSeed;
		sample = CmdArgs.INSTANCE.poly_sample_QT;
		h2 = CmdArgs.INSTANCE.simuHsq;

		E = Math.sqrt(1 - h2);

		rnd = new RandomDataImpl();
		rnd.reSeed(seed);

		freq = new double[M];
		DPrime = new double[M - 1];

		Arrays.fill(freq, CmdArgs.INSTANCE.polyFreq);
		Arrays.fill(DPrime, ld);
		LD = CalculateDprime(freq, DPrime);

		genotype = new double[sample][M];
		phenotype = new double[sample];
		BV = new double[sample];

		Logger.printUserLog("Simulation polygenic model for quantitative traits.");
		Logger.printUserLog("seed: " + seed);

		Logger.printUserLog("MAF: " + CmdArgs.INSTANCE.polyFreq);
		Logger.printUserLog("Marker: " + M);
		Logger.printUserLog("Null marker: " + M);
		if (CmdArgs.INSTANCE.polyEffectFlag)
		{
			Logger.printUserLog("genetic effect file: "
					+ CmdArgs.INSTANCE.polyEffectFile);
		} 
		else
		{
			Logger.printUserLog("Uniform Effect: " + U);
		}
		Logger.printUserLog("LD: " + ld);
		Logger.printUserLog("Sample size: " + sample);
		Logger.printUserLog("h2: " + h2);
		Logger.printUserLog("out: " + CmdArgs.INSTANCE.out);
		Logger.printUserLog("");
	}

	public static void main(String[] args)
	{
		CmdArgs.INSTANCE.parse(args);
	}

	public void generateSample()
	{
		GenerateSampleNoSelection();
		if (CmdArgs.INSTANCE.makebedFlag)
		{
			writeBFile();
		} 
		else
		{
			writeFile();
		}
	}

	public void GenerateSampleNoSelection()
	{
		DecimalFormat fmt = new DecimalFormat("#.###E0");

		int count = 0;
		RealMatrix effect = null;
		if (CmdArgs.INSTANCE.polyEffectFlag)
		{
			effect = readEffects();
		} 
		else
		{
			effect = GenerateEffects();
		}

		while (count < sample)
		{
			RealMatrix chr = SampleChromosome();
			RealMatrix res = chr.transpose().multiply(effect);

			double bv = res.getEntry(0, 0);

			BV[count] = bv;
			genotype[count] = chr.getColumn(0);

			count++;
		}

		double vg = StatUtils.variance(BV);
		double ve = vg * (1 - h2) / h2;
		double E = Math.sqrt(ve);
		Logger.printUserLog("Vg=" + fmt.format(vg));
		for (int i = 0; i < sample; i++)
		{
			phenotype[i] = BV[i] + rnd.nextGaussian(0, E);
		}
		double vp = StatUtils.variance(phenotype);
		Logger.printUserLog("Vp=" + fmt.format(vp));
		Logger.printUserLog("total individuals visited (no selection): "
				+ count + "\n");
	}

	public RealMatrix GenerateEffects()
	{

		double[] effect = new double[M];
		if (U)
		{
			for (int i = 0; i < M - M_null; i++)
			{
				double sigma_b = Math.sqrt((vy * h2)
						/ ((M - M_null) * 2 * freq[i] * (1 - freq[i])));
				effect[i] = sigma_b;
			}
		} 
		else
		{
			for (int i = 0; i < M - M_null; i++)
			{
				effect[i] = rnd.nextGaussian(0, 1);
			}
		}

		if (CmdArgs.INSTANCE.simuOrderFlag)
		{
			Arrays.sort(effect);
		}

		PrintWriter eff = null;
		try
		{
			eff = new PrintWriter(new BufferedWriter(new FileWriter(CmdArgs.INSTANCE.out
					+ ".rnd")));
		} 
		catch (IOException e)
		{
			e.printStackTrace();
		}

		for (int i = 0; i < effect.length; i++)
		{
			eff.println("rs" + i + " " + A2 + " " + effect[i]);
		}
		eff.close();

		RealMatrix Eff = new Array2DRowRealMatrix(effect);
		getMean(Eff);
		return Eff;

	}

	public RealMatrix SampleChromosome()
	{

		double[] g = new double[freq.length];
		double[][] v = new double[freq.length][2];
		for (int i = 0; i < freq.length; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				double r = rnd.nextUniform(0, 1);
				double z = r < freq[i] ? 0 : 1;
				if (i == 0)
				{
					v[i][j] = z;
				} 
				else
				{
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

	public double getMean(RealMatrix effect)
	{
		int local_sample = 10000;
		double[] P = new double[local_sample];
		int c = 0;
		while (c < local_sample)
		{

			RealMatrix chr = SampleChromosome();
			RealMatrix res = chr.transpose().multiply(effect);
			P[c++] = res.getEntry(0, 0) + rnd.nextGaussian(0, E);

		}

		Logger.printUserLog("Mean: " + StatUtils.mean(P) + ", Variance: "
				+ StatUtils.variance(P));
		return StatUtils.mean(P);
	}

	public void writeBFile()
	{
		DataOutputStream bedout = null;
		PrintWriter fam = null;
		PrintWriter bim = null;
		PrintWriter phe = null;
		PrintWriter cov = null;
		PrintWriter geno = null;
		try
		{
			bedout = new DataOutputStream(new FileOutputStream(CmdArgs.INSTANCE.out + ".bed"));

			fam = new PrintWriter(new BufferedWriter(new FileWriter(CmdArgs.INSTANCE.out
					+ ".fam")));
			bim = new PrintWriter(new BufferedWriter(new FileWriter(CmdArgs.INSTANCE.out
					+ ".bim")));

			phe = new PrintWriter(new BufferedWriter(new FileWriter(CmdArgs.INSTANCE.out
					+ ".phe")));
			cov = new PrintWriter(new BufferedWriter(new FileWriter(CmdArgs.INSTANCE.out
					+ ".cov")));
			geno = new PrintWriter(new BufferedWriter(new FileWriter(CmdArgs.INSTANCE.out
					+ ".add")));
		} 
		catch (IOException e)
		{
			e.printStackTrace();
		}

		for (int i = 0; i < genotype.length; i++)
		{
			fam.print("sample_" + i + " ");
			fam.print(1 + " ");
			fam.print(0 + " ");
			fam.print(0 + " ");
			fam.print(1 + " ");

			fam.println(phenotype[i] + " ");

			phe.println("sample_" + i + " " + 1 + " " + phenotype[i]);
			cov.println("sample_" + i + " " + 1 + " " + BV[i]);
		}

		for (int i = 0; i < genotype.length; i++)
		{
			for (int j = 0; j < genotype[i].length; j++)
			{
				geno.print(((int) genotype[i][j] + 1) + " ");
			}
			geno.println();
		}

		geno.close();

		try
		{
			bedout.writeByte(byte1);
			bedout.writeByte(byte2);
			bedout.writeByte(byte3);
			for (int i = 0; i < M; i++)
			{
				byte gbyte = 0;
				int idx = 0;
				for (int j = 0; j < sample; j++)
				{
					int g = (int) genotype[j][i] + 1;
					switch (g)
					{
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

					if (j != (sample - 1))
					{
						if (idx == 4)
						{
							bedout.writeByte(gbyte);
							gbyte = 0;
							idx = 0;
						}
					} 
					else
					{
						bedout.writeByte(gbyte);
					}
				}
			}
			bedout.close();
		} 
		catch (IOException e)
		{
			e.printStackTrace();
		}

		for (int i = 0; i < M; i++)
		{
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

	public void writeFile()
	{
		PrintWriter pedout = null;
		PrintWriter map = null;
		PrintWriter cov = null;
		PrintWriter geno = null;
		try
		{
			pedout = new PrintWriter(new BufferedWriter(new FileWriter(CmdArgs.INSTANCE.out
					+ ".ped")));
			map = new PrintWriter(new BufferedWriter(new FileWriter(CmdArgs.INSTANCE.out
					+ ".map")));
			cov = new PrintWriter(new BufferedWriter(new FileWriter(CmdArgs.INSTANCE.out
					+ ".cov")));
			geno = new PrintWriter(new BufferedWriter(new FileWriter(CmdArgs.INSTANCE.out
					+ ".add")));
		} 
		catch (IOException e)
		{
			Logger.handleException(e,
					"An exception occurred when writing files.");
		}

		for (int i = 0; i < genotype.length; i++)
		{

			pedout.print("sample_" + i + " ");
			pedout.print(1 + " ");
			pedout.print(0 + " ");
			pedout.print(0 + " ");
			pedout.print(1 + " ");
			pedout.print(1 + " ");

			for (int j = 0; j < genotype[i].length; j++)
			{
				int g = (int) genotype[i][j];
				if (g == -1)
				{
					pedout.print(A1 + " " + A1 + "  ");
				} else if (g == 0)
				{
					pedout.print(A1 + " " + A2 + "  ");
				} else
				{
					pedout.print(A2 + " " + A2 + "  ");
				}
			}
			pedout.println();

			cov.print("sample_" + i + " " + 1 + " " + 1 + " ");
			cov.print(BV[i] + " ");
			cov.println(phenotype[i]);
		}

		for (int i = 0; i < M; i++)
		{
			map.print(1 + " ");
			map.print("rs" + i + " ");
			map.print(i / (M * 1.0) + " ");
			map.println(i * 100);
		}

		for (int i = 0; i < genotype.length; i++)
		{
			for (int j = 0; j < genotype[i].length; j++)
			{
				geno.print(((int) genotype[i][j] + 1) + " ");
			}
			geno.println();
		}

		pedout.close();
		map.close();
		cov.close();
		geno.close();
	}

	public double[] CalculateDprime(double[] f, double[] dprime)
	{

		double[] D = new double[dprime.length];

		for (int i = 0; i < D.length; i++)
		{
			if (dprime[i] > 0)
			{
				D[i] = dprime[i]
						* Math.min(f[i] * (1 - f[i + 1]), f[i + 1] * (1 - f[i]));
			} 
			else
			{
				D[i] = dprime[i]
						* Math.min(f[i] * f[i + 1], (1 - f[i]) * (1 - f[i + 1]));
			}
		}

		return D;
	}

	public RealMatrix readEffects()
	{
		BufferedReader reader = FileUtil
				.FileOpen(CmdArgs.INSTANCE.polyEffectFile);
		double[] effect = new double[M];
		int c = 0;
		String line = null;
		try
		{
			while ((line = reader.readLine()) != null)
			{
				line.trim();
				String[] l = line.split(ConstValues.WHITESPACE_DELIMITER);
				if (l.length < 1)
					continue;
				if (c < (M - M_null))
				{
					effect[c++] = Double.parseDouble(l[0]);
				}
			}
		} 
		catch (IOException e)
		{
			Logger.handleException(e,
					"An exception occurred when reading the poly-effect file '"
							+ CmdArgs.INSTANCE.polyEffectFile + "'.");
		}
		RealMatrix Eff = new Array2DRowRealMatrix(effect);

		if (CmdArgs.INSTANCE.simuOrderFlag)
		{
			Arrays.sort(effect);
		}

		PrintWriter eff = null;
		try
		{
			eff = new PrintWriter(new BufferedWriter(new FileWriter(CmdArgs.INSTANCE.out
					+ ".rnd")));
		} 
		catch (IOException e)
		{
			e.printStackTrace();
		}

		for (int i = 0; i < effect.length; i++)
		{
			eff.println("rs" + i + " " + A2 + " " + effect[i]);
		}
		eff.close();

		return Eff;
	}
}
