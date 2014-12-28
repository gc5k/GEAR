package gear.simulation;

import gear.CmdArgs;
import gear.ConstValues;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.pop.PopStat;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
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

public class SimuPolyCC
{

	private RandomDataImpl rnd;
	private long seed = 2011;

	RealMatrix effect = null;

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

	private double h2 = 0.5;
	private double E = 0.5;
	private double K = 0.05;

	private String A1 = "A";
	private String A2 = "C";
	public static StringBuilder LOG = new StringBuilder();

	public SimuPolyCC()
	{

		M = CmdArgs.INSTANCE.polyLoci;
		M_null = CmdArgs.INSTANCE.polyLociNull;
		U = CmdArgs.INSTANCE.polyU;
		ld = CmdArgs.INSTANCE.polyLD;
		seed = CmdArgs.INSTANCE.simuSeed;
		sample = CmdArgs.INSTANCE.simuCC[0] + CmdArgs.INSTANCE.simuCC[1];
		N_case = CmdArgs.INSTANCE.simuCC[0];
		N_control = CmdArgs.INSTANCE.simuCC[1];
		h2 = CmdArgs.INSTANCE.simuHsq;
		K = CmdArgs.INSTANCE.simuK;

		if (Math.abs(h2-1) < 0.001) 
		{
			Logger.printUserLog("'h2 = " + h2 +"' is too big.  Now it has been reset to " + (h2 * 0.99));
			h2 *= 0.99;
		}

		E = Math.sqrt(1 - h2);// will be rescaled in scaling();

		rnd = new RandomDataImpl();
		rnd.reSeed(seed);

		freq = new double[M];
		DPrime = new double[M - 1];

		Arrays.fill(freq, CmdArgs.INSTANCE.polyFreq);
		Arrays.fill(DPrime, ld);
		LD = PopStat.CalcLDfromDPrime(freq, DPrime);

		genotype = new double[N_case + N_control][M];
		phenotype = new double[N_case + N_control];
		BV = new double[N_case + N_control];
		risk = new double[N_case + N_control];

		Calendar calendar = Calendar.getInstance();
		LOG.append("\nThe analysis was implemented at: " + calendar.getTime()
				+ "\n");
		LOG.append("seed: " + seed + "\n");
		LOG.append("MAF: " + CmdArgs.INSTANCE.polyFreq + "\n");
		LOG.append("Marker: " + M + "\n");
		LOG.append("Null Marker: " + M_null + "\n");
		if (CmdArgs.INSTANCE.polyEffectFlag)
		{
			effect = readEffects();
			LOG.append("Genetic effect file: "
					+ CmdArgs.INSTANCE.polyEffectFile + "\n");
		} 
		else
		{
			effect = generateEffects();
			LOG.append("Uniform Effect: " + U + "\n");
		}
		LOG.append("LD (Lewontin's Dprime): " + ld + "\n"); 
		LOG.append("Sample size: " + sample + "\n");
		LOG.append("Case: " + N_case + "\n");
		LOG.append("Control: " + N_control + "\n");
		LOG.append("K (prevalence): " + K + "\n");
		LOG.append("h2: " + h2 + "\n");
		LOG.append("out: " + CmdArgs.INSTANCE.out + "\n");
		if (CmdArgs.INSTANCE.makebedFlag)
		{
			LOG.append("Generated files in bed format.");
		}
		LOG.append("\n");
		
		scaling();
		generateSampleSelection();
		if (CmdArgs.INSTANCE.makebedFlag)
		{
			writeBFile();
		}
		else
		{
			writeFile();
		}
		writeLog();
	}

	private void scaling()
	{
		int TestRun = N_case+N_control;
		double[] bv = new double[TestRun];
		Logger.printUserLog("Rescaling with " + TestRun + " runs.");

		for(int i = 0; i < TestRun; i++)
		{
			RealMatrix chr = SampleChromosome();
			RealMatrix genoEff = chr.transpose().multiply(effect);	
			bv[i] = genoEff.getEntry(0, 0);
		}

		double vg = StatUtils.variance(bv);
		double m = StatUtils.mean(bv);
		double ve = vg*(1-h2)/h2;
		E = Math.sqrt(ve);

		Logger.printUserLog("Mean of breeding values: " + m);
		Logger.printUserLog("Vg: "+ vg);
		Logger.printUserLog("Rescale residual to " + E);
		Logger.printUserLog("Vp: " + vg/h2);

		double sd = Math.sqrt(E*E/(1-h2));
		NormalDistributionImpl norm = new NormalDistributionImpl(0, sd);

		double cnt=0;
		for(int i = 0; i < TestRun; i++)
		{
			double L = bv[i] + rnd.nextGaussian(0, E);
			double liability = 0;
			try
			{
				liability = norm.cumulativeProbability(L);
			} 
			catch (MathException e)
			{
				e.printStackTrace();
			}
			
			if ( (1 - liability) < K)
			{
				cnt++;
			}
		}
		Logger.printUserLog("Emperical K: " + cnt/(TestRun));
	}

	public void generateSampleSelection()
	{
		int count_case = 0;
		int count_control = 0;
		int count = 0;

		double sd = Math.sqrt(E*E/(1-h2));
		NormalDistributionImpl norm = new NormalDistributionImpl(0, sd);

		int cnt = 0;
		while (count < sample)
		{
			cnt++;
			RealMatrix chr = SampleChromosome();
			RealMatrix genoEff = chr.transpose().multiply(effect);

			double bv = genoEff.getEntry(0, 0);
			double L = bv + rnd.nextGaussian(0, E);
			double liability = 0;

			try
			{
				liability = norm.cumulativeProbability(L);
			} 
			catch (MathException e)
			{
				e.printStackTrace();
			}

			if (1 - liability < K)
			{
				if (count_case < N_case)
				{
					BV[count_case] = bv;
					phenotype[count_case] = L;
					genotype[count_case] = chr.getColumn(0);
					risk[count_case] = liability;
					count_case++;
				} 
				else
				{
					continue;
				}
			} 
			else
			{
				if (count_control < N_control)
				{
					BV[N_case + count_control] = bv;
					phenotype[N_case + count_control] = L;
					genotype[N_case + count_control] = chr.getColumn(0);
					risk[N_case + count_control] = liability;
					count_control++;
				} 
				else
				{
					continue;
				}
			}
			count++;
		}

		double EpCnt = N_case/K;
		LOG.append("Total individuals have been visited (expected): " + cnt + " (" + EpCnt+ ")\n");
		check();
	}
	
	private void check()
	{
		double m_cs, m_ctrl, v_cs, v_cs_a, v_ctrl, v_ctrl_a;
		double[] cs = new double[N_case];
		double[] cs_a = new double[N_case];
		System.arraycopy(phenotype, 0, cs, 0, N_case);
		System.arraycopy(BV, 0, cs_a, 0, N_case);
		m_cs = StatUtils.mean(cs);
		v_cs = StatUtils.variance(cs);
		v_cs_a = StatUtils.variance(cs_a);

		double[] ctrl = new double[N_control];
		double[] ctrl_a = new double[N_control];
		System.arraycopy(phenotype, N_case, ctrl, 0, N_control);
		System.arraycopy(BV, N_case, ctrl_a, 0, N_control);
		m_ctrl = StatUtils.mean(ctrl);
		v_ctrl = StatUtils.variance(ctrl);
		v_ctrl_a = StatUtils.variance(ctrl_a);

		double sd_p = Math.sqrt(E*E/(1-h2));
		NormalDistributionImpl norm = new NormalDistributionImpl(0, 1);
		double threshold = 0;
		double x = 0; 
		try
		{
			x = norm.inverseCumulativeProbability(1-K);
			threshold = norm.density(x);
		}
		catch (MathException e)
		{
			e.printStackTrace();
		}
		double selectionIntensity = threshold/K;
		LOG.append("Selection intensity: " + selectionIntensity + "\n");

		double Em_cs, Em_ctrl, Ev_cs, Ev_cs_a, Ev_ctrl, Ev_ctrl_a;
		Em_cs = selectionIntensity * sd_p;
		Em_ctrl = -1* Em_cs * K/(1-K);

		double threshold_ctrl = 0;
		double x_ctrl = 0; 
		try
		{
			x_ctrl = norm.inverseCumulativeProbability(K);
			threshold_ctrl = norm.density(x_ctrl);
		}
		catch (MathException e)
		{
			e.printStackTrace();
		}

		double selectionIntensity_ctrl = threshold_ctrl/(1-K);
		Ev_cs = selectionIntensity_ctrl * (selectionIntensity_ctrl - x_ctrl) * sd_p * sd_p;
		Ev_ctrl = selectionIntensity * (selectionIntensity - x) * sd_p * sd_p;

		Ev_cs_a = (1 - h2*selectionIntensity * (selectionIntensity - x)) * sd_p * sd_p * h2;
		Ev_ctrl_a = (1 - h2*selectionIntensity_ctrl * (selectionIntensity_ctrl - x_ctrl)) * sd_p * sd_p * h2;
		
		LOG.append("Expected mean of control phenotype (observed): " + Em_ctrl + " (" + m_ctrl + ")\n");
		LOG.append("Expected variance of controls' phenotype (observed): " + Ev_ctrl + " (" + v_ctrl + ")\n");
		LOG.append("Expected variance of controls' liability (observed): " + Ev_ctrl_a + " (" + v_ctrl_a + ")\n");

		LOG.append("Expected mean of case phenotype (observed): " + Em_cs + " (" + m_cs + ")\n");
		LOG.append("Expected variance of cases' phenotype (observed): " + Ev_cs + " (" + v_cs + ")\n");
		LOG.append("Expected variance of cases' liability (observed): " + Ev_cs_a + " (" + v_cs_a + ")\n");

	}

	public static void main(String[] args)
	{
		CmdArgs.INSTANCE.parse(args);
	}

	public void writeLog()
	{
		Calendar calendar = Calendar.getInstance();
		LOG.append("\nThe analysis was completed at: " + calendar.getTime()
				+ "\n");
		PrintWriter log = null;
		try
		{
			log = new PrintWriter(new BufferedWriter(new FileWriter(CmdArgs.INSTANCE.out
					+ ".plog")));
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
		log.println(LOG.toString());
		log.close();

		Logger.printUserLog(LOG.toString());
	}

	private RealMatrix generateEffects()
	{

		double[] effect = new double[M];
		if (U)
		{
			for (int i = 0; i < effect.length - M_null; i++)
			{
				double sigma_b = Math.sqrt((h2)
						/ ((M - M_null) * 2 * freq[i] * (1 - freq[i])));
				effect[i] = sigma_b;
			}
		}
		else
		{
			for (int i = 0; i < effect.length - M_null; i++)
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
			if (i < N_case)
			{
				fam.print("case_" + i + " ");
				fam.print(1 + " ");
				fam.print(0 + " ");
				fam.print(0 + " ");
				fam.print(1 + " ");

				fam.println(2 + " ");
			} 
			else
			{
				fam.print("control_" + i + " ");
				fam.print(1 + " ");
				fam.print(0 + " ");
				fam.print(0 + " ");
				fam.print(1 + " ");

				fam.println(1 + " ");
			}

			if (i < N_case)
			{
				phe.println("case_" + i + " " + 1 + " " + 2 + " ");
				cov.print("case_" + i + " " + 1 + " " + 2 + " ");
			} 
			else
			{
				phe.println("control_" + i + " " + 1 + " " + 1 + " ");
				cov.print("control_" + i + " " + 1 + " " + 1 + " ");
			}

			cov.print(risk[i] + " ");
			cov.println(BV[i] + " ");
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
			bedout.writeByte(ConstValues.PLINK_BED_BYTE1);
			bedout.writeByte(ConstValues.PLINK_BED_BYTE2);
			bedout.writeByte(ConstValues.PLINK_BED_BYTE3);
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
		PrintWriter phe = null;
		PrintWriter cov = null;
		PrintWriter geno = null;

		try
		{
			pedout = new PrintWriter(new BufferedWriter(new FileWriter(CmdArgs.INSTANCE.out
					+ ".ped")));
			map = new PrintWriter(new BufferedWriter(new FileWriter(CmdArgs.INSTANCE.out
					+ ".map")));
			phe = new PrintWriter(new BufferedWriter(new FileWriter(CmdArgs.INSTANCE.out
					+ ".phe")));
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
			if (i < N_case)
			{
				pedout.print("case_" + i + " ");
				pedout.print(1 + " ");
				pedout.print(0 + " ");
				pedout.print(0 + " ");
				pedout.print(1 + " ");

				pedout.print(2 + " ");
			} 
			else
			{
				pedout.print("control_" + i + " ");
				pedout.print(1 + " ");
				pedout.print(0 + " ");
				pedout.print(0 + " ");
				pedout.print(1 + " ");

				pedout.print(1 + " ");
			}

			for (int j = 0; j < genotype[i].length; j++)
			{
				int g = (int) genotype[i][j];
				if (g == -1)
				{
					pedout.print(A1 + " " + A1 + "  ");
				} 
				else if (g == 0)
				{
					pedout.print(A1 + " " + A2 + "  ");
				} 
				else
				{
					pedout.print(A2 + " " + A2 + "  ");
				}
			}
			pedout.println();

			if (i < N_case)
			{
				phe.println("case_" + i + " " + 1 + " " + 2 + " ");
				cov.print("case_" + i + " " + 1 + " " + 2 + " ");
			} 
			else
			{
				phe.println("control_" + i + " " + 1 + " " + 1 + " ");
				cov.print("control_" + i + " " + 1 + " " + 1 + " ");
			}

			cov.print(risk[i] + " ");
			cov.println(BV[i] + " ");
		}

		for (int i = 0; i < genotype.length; i++)
		{
			for (int j = 0; j < genotype[i].length; j++)
			{
				geno.print(((int) genotype[i][j] + 1) + " ");
			}
			geno.println();
		}

		for (int i = 0; i < M; i++)
		{
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
