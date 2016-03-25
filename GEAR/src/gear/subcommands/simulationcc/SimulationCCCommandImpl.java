package gear.subcommands.simulationcc;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistributionImpl;
import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.random.RandomDataImpl;
import org.apache.commons.math.stat.StatUtils;

import gear.ConstValues;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.Sample;
import gear.util.pop.PopStat;

public class SimulationCCCommandImpl extends CommandImpl
{

	@Override
	public void execute(CommandArguments cmdArgs)
	{
		ccArgs = (SimulationCCCommandArguments)cmdArgs;

		sample_cs = ccArgs.getCaseSize();
		sample_ctrl = ccArgs.getControlSize();

		h2 = ccArgs.getHsq();
		M = ccArgs.getMarkerNum();
		nullM = ccArgs.getNullMarkerNum();
		seed = ccArgs.getSeed();
		rnd.reSeed(seed);
		risk = new double[sample_cs + sample_ctrl];
		genotype = new double[sample_cs + sample_ctrl][M];
		phenotype = new double[sample_cs + sample_ctrl];
		BV = new double[sample_cs + sample_ctrl];

		getFreq();
		getEffect();
		getDPrime();
		calLD();
		scaling();
		generateSampleSelection();
		
		if (ccArgs.isMakeBed())
		{
			writeBFile();
		} 
		else
		{
			writeFile();
		}

	}

	private void scaling()
	{
		int TestRun = sample_cs+sample_ctrl;
		double[] bv = new double[TestRun];
		Logger.printUserLog("Rescaling with " + TestRun + " runs.");

		RealMatrix Meffect = new Array2DRowRealMatrix(effect);

		for(int i = 0; i < TestRun; i++)
		{
			RealMatrix chr = SampleChromosome();
			RealMatrix genoEff = chr.transpose().multiply(Meffect);	
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
			
			if ( (1 - liability) < ccArgs.getK())
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
		RealMatrix Meffect = new Array2DRowRealMatrix(effect);

		double sd = Math.sqrt(E*E/(1-h2));
		NormalDistributionImpl norm = new NormalDistributionImpl(0, sd);

		int cnt = 0;
		while (count < (sample_cs + sample_ctrl))
		{
			cnt++;
			RealMatrix chr = SampleChromosome();
			RealMatrix genoEff = chr.transpose().multiply(Meffect);

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

			if (1 - liability < ccArgs.getK())
			{
				if (count_case < sample_cs)
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
				if (count_control < sample_ctrl)
				{
					BV[sample_cs + count_control] = bv;
					phenotype[sample_cs + count_control] = L;
					genotype[sample_cs + count_control] = chr.getColumn(0);
					risk[sample_cs + count_control] = liability;
					count_control++;
				} 
				else
				{
					continue;
				}
			}
			count++;
		}

		double EpCnt = sample_cs/ccArgs.getK();
		Logger.printUserLog("Total individuals have been visited (expected): " + cnt + " (" + EpCnt+ ")");
		check();
	}

	private void check()
	{
		double m_cs, m_ctrl, v_cs, v_cs_a, v_ctrl, v_ctrl_a;
		double[] cs = new double[sample_cs];
		double[] cs_a = new double[sample_cs];
		System.arraycopy(phenotype, 0, cs, 0, sample_cs);
		System.arraycopy(BV, 0, cs_a, 0, sample_cs);
		m_cs = StatUtils.mean(cs);
		v_cs = StatUtils.variance(cs);
		v_cs_a = StatUtils.variance(cs_a);

		double[] ctrl = new double[sample_ctrl];
		double[] ctrl_a = new double[sample_ctrl];
		System.arraycopy(phenotype, sample_cs, ctrl, 0, sample_ctrl);
		System.arraycopy(BV, sample_cs, ctrl_a, 0, sample_ctrl);
		m_ctrl = StatUtils.mean(ctrl);
		v_ctrl = StatUtils.variance(ctrl);
		v_ctrl_a = StatUtils.variance(ctrl_a);

		double sd_p = Math.sqrt(E*E/(1-h2));
		NormalDistributionImpl norm = new NormalDistributionImpl(0, 1);
		double threshold = 0;
		double x = 0; 
		try
		{
			x = norm.inverseCumulativeProbability(1-ccArgs.getK());
			threshold = norm.density(x);
		}
		catch (MathException e)
		{
			e.printStackTrace();
		}
		double selectionIntensity = threshold/ccArgs.getK();
		Logger.printUserLog("Selection intensity: " + selectionIntensity);

		double Em_cs, Em_ctrl, Ev_cs, Ev_cs_a, Ev_ctrl, Ev_ctrl_a;
		Em_cs = selectionIntensity * sd_p;
		Em_ctrl = -1* Em_cs * ccArgs.getK()/(1-ccArgs.getK());

		double threshold_ctrl = 0;
		double x_ctrl = 0; 
		try
		{
			x_ctrl = norm.inverseCumulativeProbability(ccArgs.getK());
			threshold_ctrl = norm.density(x_ctrl);
		}
		catch (MathException e)
		{
			e.printStackTrace();
		}

		double selectionIntensity_ctrl = threshold_ctrl/(1-ccArgs.getK());
		Ev_cs = selectionIntensity_ctrl * (selectionIntensity_ctrl - x_ctrl) * sd_p * sd_p;
		Ev_ctrl = selectionIntensity * (selectionIntensity - x) * sd_p * sd_p;

		Ev_cs_a = (1 - h2*selectionIntensity * (selectionIntensity - x)) * sd_p * sd_p * h2;
		Ev_ctrl_a = (1 - h2*selectionIntensity_ctrl * (selectionIntensity_ctrl - x_ctrl)) * sd_p * sd_p * h2;
	
		Logger.printUserLog("Observed mean of control phenotype (expected): " + m_ctrl + " (" + Em_ctrl + ")");
		Logger.printUserLog("Observed variance of controls' phenotype (expected): " + v_ctrl + " (" + Ev_ctrl + ")");
		Logger.printUserLog("Observed variance of controls' liability (expected): " + v_ctrl_a + " (" + Ev_ctrl_a + ")");

		Logger.printUserLog("Observed mean of case phenotype (expected): " + m_cs + " (" + Em_cs + ")");
		Logger.printUserLog("Observed variance of cases' phenotype (expected): " + v_cs + " (" + Ev_cs + ")");
		Logger.printUserLog("Observed variance of cases' liability (expected): " + v_cs_a + " (" + Ev_cs_a + ")");
	}

	public RealMatrix SampleChromosome()
	{
		double[] g = new double[M];
		double[][] v = new double[M][2];
		for (int i = 0; i < M; i++)
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

	
	private void getFreq()
	{
		freq = new double[M];

		if (ccArgs.isPlainFreq())
		{
			Arrays.fill(freq, ccArgs.getFreq());
		}
		else if (ccArgs.isUnifFreq())
		{
			for(int i = 0; i < M; i++)
			{
				freq[i] = rnd.nextUniform(ccArgs.getFreqRangeLow(), ccArgs.getFreqRangeHigh());
			}
		}
		else if (ccArgs.isFreqFile())
		{
			BufferedReader reader = FileUtil.FileOpen(ccArgs.getFreqFile());
			int c = 0;
			String line = null;
			try
			{
				while ((line = reader.readLine()) != null)
				{
					if(c >= M)
					{
						Logger.printUserLog("Have already read " + M + " allelic frequencies.  Ignore the rest of the content in '" + ccArgs.getFreqFile() + "'.");
					}
					line.trim();
					String[] l = line.split(ConstValues.WHITESPACE_DELIMITER);
					if (l.length < 1) continue;
					freq[c++] = Double.parseDouble(l[0]);
				}
				reader.close();
			}
			catch (IOException e)
			{
				Logger.handleException(e,
						"An exception occurred when reading the frequency file '"
								+ ccArgs.getFreqFile() + "'.");
			}
		}
	}

	private void getEffect()
	{
		effect = new double[M];
		
		Sample.setSeed(seed);
		
		int[] idx = Sample.SampleIndex(0, M-1, M-nullM);
		Arrays.sort(idx);

		if(ccArgs.isPlainEffect())
		{
			for(int i = 0; i < idx.length; i++) effect[idx[i]] = ccArgs.getPolyEffect();
		}
		else if (ccArgs.isPolyEffect())
		{
			for (int i = 0; i < idx.length; i++)
			{
				effect[idx[i]] = rnd.nextGaussian(0, 1);
			}
		}
		else if (ccArgs.isPolyEffectSort())
		{
			NormalDistributionImpl ndImpl = new NormalDistributionImpl();
			ndImpl.reseedRandomGenerator(ccArgs.getSeed());
			for (int i = 0; i < idx.length; i++)
			{
				try
				{
					effect[idx[i]] = ndImpl.inverseCumulativeProbability((i+0.5)/(idx.length));
				}
				catch (MathException e)
				{
					e.printStackTrace();
				}
			}
		}
		else if (ccArgs.isPolyEffectFile())
		{
			BufferedReader reader = FileUtil.FileOpen(ccArgs.getPolyEffectFile());
			int c = 0;
			String line = null;
			try
			{
				while ((line = reader.readLine()) != null)
				{
					if(c >= M)
					{
						Logger.printUserLog("Have already read " + M + " allelic effects.  Ignore the rest of the content in '" + ccArgs.getPolyEffectFile() + "'.");
						break;
					}

					line.trim();
					String[] l = line.split(ConstValues.WHITESPACE_DELIMITER);
					if (l.length < 1) continue;
					effect[c++] = Double.parseDouble(l[0]);
				}
				reader.close();
			}
			catch (IOException e)
			{
				Logger.handleException(e,
						"An exception occurred when reading the frequency file '"
								+ ccArgs.getPolyEffectFile() + "'.");
			}
		}
	}

	private void getDPrime()
	{
		dprime = new double[M-1];
		if(ccArgs.isPlainLD())
		{
			Arrays.fill(dprime, ccArgs.getLD());
		}
		else if (ccArgs.isRandLD())
		{
			for(int i = 0; i < dprime.length; i++)
			{
				dprime[i] = rnd.nextUniform(ccArgs.getFreqRangeLow(), ccArgs.getFreqRangeHigh());
			}
		}
	}

	public void calLD()
	{
		LD = PopStat.CalcLDfromDPrime(freq, dprime);
	}

	public void writeBFile()
	{
		DataOutputStream bedout = null;
		PrintWriter fam = null;
		PrintWriter bim = null;
		PrintWriter phe = null;
		PrintWriter geno = null;
		PrintWriter eff = null;

		try
		{
			bedout = new DataOutputStream(new FileOutputStream(ccArgs.getOutRoot() + ".bed"));

			fam = new PrintWriter(new BufferedWriter(new FileWriter(ccArgs.getOutRoot()
					+ ".fam")));
			bim = new PrintWriter(new BufferedWriter(new FileWriter(ccArgs.getOutRoot()
					+ ".bim")));

			phe = new PrintWriter(new BufferedWriter(new FileWriter(ccArgs.getOutRoot()
					+ ".phe")));
			geno = new PrintWriter(new BufferedWriter(new FileWriter(ccArgs.getOutRoot()
					+ ".add")));
			eff = new PrintWriter(new BufferedWriter(new FileWriter(ccArgs.getOutRoot() + ".rnd")));

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

			phe.print("sample_" + i + " " + 1 + " " + BV[i]);
			phe.print(" " + phenotype[i] + " ");

			if (i < sample_cs)
			{
				fam.println("2");
				phe.println("2");
			}
			else
			{
				fam.println("1");
				phe.println("1");
			}
		}

		for (int i = 0; i < genotype.length; i++)
		{
			for (int j = 0; j < genotype[i].length; j++)
			{
				geno.print(((int) genotype[i][j] + 1) + " ");
			}
			geno.println();
		}

		try
		{
			bedout.writeByte(ConstValues.PLINK_BED_BYTE1);
			bedout.writeByte(ConstValues.PLINK_BED_BYTE2);
			bedout.writeByte(ConstValues.PLINK_BED_BYTE3);
			for (int i = 0; i < M; i++)
			{
				byte gbyte = 0;
				int idx = 0;
				for (int j = 0; j < (sample_cs + sample_ctrl); j++)
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

					if (j != ( (sample_cs +sample_ctrl) - 1))
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
			eff.println("rs" + i + " " + A1 + " " + effect[i]);
		}

		geno.close();
		phe.close();
		bim.close();
		fam.close();
		eff.close();
	}

	public void writeFile()
	{
		PrintWriter pedout = null;
		PrintWriter map = null;
		PrintWriter phe = null;
		PrintWriter geno = null;
		PrintWriter eff = null;
		try
		{
			pedout = new PrintWriter(new BufferedWriter(new FileWriter(ccArgs.getOutRoot()
					+ ".ped")));
			map = new PrintWriter(new BufferedWriter(new FileWriter(ccArgs.getOutRoot()
					+ ".map")));
			phe = new PrintWriter(new BufferedWriter(new FileWriter(ccArgs.getOutRoot()
					+ ".phe")));
			geno = new PrintWriter(new BufferedWriter(new FileWriter(ccArgs.getOutRoot()
					+ ".add")));
			eff = new PrintWriter(new BufferedWriter(new FileWriter(ccArgs.getOutRoot() + ".rnd")));
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
			if(i < sample_cs)
			{
				pedout.print(2 + " ");
			}
			else
			{
				pedout.print(1 + " ");
			}
			pedout.print(phenotype[i] + " ");

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

			phe.print("sample_" + i + " " + 1 + " " + BV[i]);
			phe.print(" " + phenotype[i]);
			if(i < sample_cs)
			{
				phe.println(" " + 2);
			}
			else
			{
				phe.println(" " + 1);
			}
		}

		for (int i = 0; i < M; i++)
		{
			map.print(1 + " ");
			map.print("rs" + i + " ");
			map.print(i / (M * 1.0) + " ");
			map.println(i * 100);
			eff.println("rs" + i + " " + A1 + " " + effect[i]);
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
		phe.close();
		geno.close();
		eff.close();
	}
	

	private SimulationCCCommandArguments ccArgs;
	
	private RandomDataImpl rnd = new RandomDataImpl();
	private long seed;

	private int M;
	private int nullM;
	private int sample_cs;
	private int sample_ctrl;

	private double E;
	private double[][] genotype;
	private double[] BV;
	private double[] phenotype;

	private double[] effect;
	private double[] freq;

	private double[] dprime;
	private double[] LD;
	private double[] risk;

	private double h2;

	private String A1 = "A";
	private String A2 = "C";

}
