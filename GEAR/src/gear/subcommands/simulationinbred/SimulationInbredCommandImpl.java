package gear.subcommands.simulationinbred;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
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

public class SimulationInbredCommandImpl extends CommandImpl
{

	@Override
	public void execute(CommandArguments cmdArgs)
	{
		inbredArgs = (SimulationInbredCommandArguments)cmdArgs;

		sample = inbredArgs.getSampleSize();
		h2 = inbredArgs.getHsq();
		M = inbredArgs.getMarkerNum();
		nullM = inbredArgs.getNullMarkerNum();
		seed = inbredArgs.getSeed();
		rnd.reSeed(seed);
		rep = inbredArgs.getRep();

		getFreq();
		getEffect();
		getDPrime();
		calLD();
		generateSampleNoSelection();

		if (inbredArgs.isMakeBed())
		{
			writeBFile();
		}
		else
		{
			writeFile();
		}
		writeEffFile();

	}

	private void generateSampleNoSelection()
	{
		DecimalFormat fmt = new DecimalFormat("#.###E0");

		RealMatrix Meffect = new Array2DRowRealMatrix(effect);
		genotype = new double[sample][M];
		phenotype = new double[sample][rep];
		BV = new double[sample];

		for (int i = 0; i < sample; i++)
		{
			RealMatrix chr = SampleChromosome();
			RealMatrix genoEff = chr.transpose().multiply(Meffect);

			double bv = genoEff.getEntry(0, 0);
			BV[i] = bv;
			genotype[i] = chr.getColumn(0);
		}

		double H2=inbredArgs.getHsq();
		if (H2 == 0)
		{
			Arrays.fill(BV, 0);
		}

		double vg = StatUtils.variance(BV);
		//rescale the phenotype to get the heritability and residual
		double ve = H2 == 0 ? 1:vg * (1 - H2) / H2;
		double E = Math.sqrt(ve);
		Logger.printUserLog("Vg=" + fmt.format(vg));
		for (int i = 0; i < rep; i++)
		{
			double[] pv=new double[sample];
			for (int j = 0; j < sample; j++)
			{
				phenotype[j][i] = BV[j] + rnd.nextGaussian(0, E);
				pv[j] = phenotype[j][i];
			}
			double Vp=StatUtils.variance(pv);
			Logger.printUserLog("Vp=" + fmt.format(Vp) + "; hsq=" + fmt.format(vg/Vp) + " for replicate " + (i+1));
		}
		Logger.printUserLog("Total individuals visited (no selection): "
				+ BV.length + "\n");
	}

	public RealMatrix SampleChromosome()
	{
		double[] g = new double[M];
		double[][] v = new double[M][2];
		for (int i = 0; i < M; i++)
		{
			for (int j = 0; j < 1; j++)
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
				v[i][1]=v[i][0];
			}
			g[i] = v[i][0] + v[i][1] - 1;
		}

		RealMatrix chr = new Array2DRowRealMatrix(g);
		return chr;
	}

	private void getFreq()
	{
		freq = new double[M];

		if (inbredArgs.isPlainFreq())
		{
			Arrays.fill(freq, inbredArgs.getFreq());
		}
		else if (inbredArgs.isUnifFreq())
		{
			for(int i = 0; i < M; i++)
			{
				freq[i] = rnd.nextUniform(inbredArgs.getFreqRangeLow(), inbredArgs.getFreqRangeHigh());
			}
		}
		else if (inbredArgs.isFreqFile())
		{
			BufferedReader reader = FileUtil.FileOpen(inbredArgs.getFreqFile());
			int c = 0;
			String line = null;
			try
			{
				while ((line = reader.readLine()) != null)
				{
					if(c >= M)
					{
						Logger.printUserLog("Have already read " + M + " allelic frequencies. Ignore the rest of the content in '" + inbredArgs.getFreqFile() + "'.");
						break;
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
								+ inbredArgs.getFreqFile() + "'.");
			}
		}
	}

	private void getEffect()
	{
		effect = new double[M];
		Sample.setSeed(seed);
		int[] idx = Sample.SampleIndex(0, M-1, M-nullM);
		Arrays.sort(idx);

		if(inbredArgs.isPlainEffect())
		{
			for(int i = 0; i < idx.length; i++) effect[idx[i]] = inbredArgs.getPolyEffect();
		}
		else if (inbredArgs.isPolyEffect())
		{
			for (int i = 0; i < idx.length; i++)
			{
				effect[idx[i]] = rnd.nextGaussian(0, 1);
			}
		}
		else if (inbredArgs.isPolyEffectSort())
		{
			NormalDistributionImpl ndImpl = new NormalDistributionImpl();
			ndImpl.reseedRandomGenerator(inbredArgs.getSeed());
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
		else if (inbredArgs.isPolyEffectFile())
		{
			BufferedReader reader = FileUtil.FileOpen(inbredArgs.getPolyEffectFile());
			int c = 0;
			String line = null;
			try
			{
				while ((line = reader.readLine()) != null)
				{
					if(c >= M)
					{
						Logger.printUserLog("Have already read " + M + " allelic effects. Ignore the rest of the content in '" + inbredArgs.getPolyEffectFile() + "'.");
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
								+ inbredArgs.getPolyEffectFile() + "'.");
			}
		}
		
		if (h2 == 0)
		{
			Arrays.fill(effect, 0);
		}
	}


	private void getDPrime()
	{
		dprime = new double[M-1];
		if(inbredArgs.isPlainLD())
		{
			Arrays.fill(dprime, inbredArgs.getLD());
		}
		else if (inbredArgs.isRandLD())
		{
			for(int i = 0; i < dprime.length; i++)
			{
				dprime[i] = rnd.nextUniform(-1, 1);
			}
		}
		else if (inbredArgs.isLDFile())
		{
			BufferedReader reader = FileUtil.FileOpen(inbredArgs.getLDFile());
			int c = 0;
			String line = null;
			try
			{
				while ((line = reader.readLine()) != null)
				{
					if (c >= dprime.length)
					{
						Logger.printUserLog("Have already read " + dprime.length + " LD. Ignore the rest of the content in '" + inbredArgs.getLDFile() + "'.");
						break;
					}

					line.trim();
					String[] l = line.split(ConstValues.WHITESPACE_DELIMITER);
					if (l.length < 1) continue;
					dprime[c++] = Double.parseDouble(l[0]);
				}
				reader.close();
			}
			catch (IOException e)
			{
				Logger.handleException(e,
						"An exception occurred when reading the ld file '"
								+ inbredArgs.getLDFile() + "'.");
			}
		}
	}

	public void calLD()
	{
		LD = PopStat.CalcLDfromDPrime(freq, dprime);
//		LD = new double[M-1];
//
//		for (int i = 0; i < LD.length; i++)
//		{
//			if (dprime[i] > 0)
//			{
//				LD[i] = dprime[i]
//						* Math.min(freq[i] * (1 - freq[i + 1]), freq[i + 1] * (1 - freq[i]));
//			} 
//			else
//			{
//				LD[i] = dprime[i]
//						* Math.min(freq[i] * freq[i + 1], (1 - freq[i]) * (1 - freq[i + 1]));
//			}
//		}
	}

	public void writeBFile()
	{
		DataOutputStream bedout = null;
		PrintWriter fam = null;
		PrintWriter bim = null;
		PrintWriter phe = null;
		PrintWriter geno = null;
		PrintWriter breed = null;

		try
		{
			bedout = new DataOutputStream(new FileOutputStream(inbredArgs.getOutRoot() + ".bed"));

			fam = new PrintWriter(new BufferedWriter(new FileWriter(inbredArgs.getOutRoot() + ".fam")));
			bim = new PrintWriter(new BufferedWriter(new FileWriter(inbredArgs.getOutRoot() + ".bim")));
			phe = new PrintWriter(new BufferedWriter(new FileWriter(inbredArgs.getOutRoot() + ".phe")));
			geno = new PrintWriter(new BufferedWriter(new FileWriter(inbredArgs.getOutRoot() + ".add")));
			breed = new PrintWriter(new BufferedWriter(new FileWriter(inbredArgs.getOutRoot() + ".breed")));
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}

		for (int i = 0; i < genotype.length; i++)
		{
			if (inbredArgs.isFID())
			{
				fam.print(inbredArgs.getFamIDPrefix() + "sample_" + i + " ");
			}
			else
			{
				fam.print("sample_" + i + " ");				
			}
			fam.print(1 + " ");
			fam.print(0 + " ");
			fam.print(0 + " ");
			fam.print(1 + " ");

			fam.println(phenotype[i][0] + " ");

			if (inbredArgs.isFID())
			{
				phe.print("sample_" + i + " " + 1);
				breed.println("sample_" + i + " " + 1 + " " + BV[i]);				
			}
			else
			{
				phe.print(inbredArgs.getFamIDPrefix() + "sample_" + i + " " + 1);
				breed.println(inbredArgs.getFamIDPrefix() + "sample_" + i + " " + 1 + " " + BV[i]);				
			}

			for (int j = 0; j < rep; j++)
			{
				phe.print(" " + phenotype[i][j]);
			}
			phe.println();
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

		geno.close();
		phe.close();
		bim.close();
		fam.close();
		breed.close();
	}

	public void writeFile()
	{
		PrintWriter pedout = null;
		PrintWriter map = null;
		PrintWriter phe = null;
		PrintWriter geno = null;
		PrintWriter breed = null;
		try
		{
			pedout = new PrintWriter(new BufferedWriter(new FileWriter(inbredArgs.getOutRoot() + ".ped")));
			map = new PrintWriter(new BufferedWriter(new FileWriter(inbredArgs.getOutRoot() + ".map")));
			phe = new PrintWriter(new BufferedWriter(new FileWriter(inbredArgs.getOutRoot() + ".phe")));
			geno = new PrintWriter(new BufferedWriter(new FileWriter(inbredArgs.getOutRoot() + ".add")));
			breed = new PrintWriter(new BufferedWriter(new FileWriter(inbredArgs.getOutRoot() + ".breed")));
		}
		catch (IOException e)
		{
			Logger.handleException(e,
					"An exception occurred when writing files.");
		}

		for (int i = 0; i < genotype.length; i++)
		{
			if (inbredArgs.isFID())
			{
				pedout.print(inbredArgs.getFamIDPrefix() + "sample_" + i + " ");				
			}
			else 
			{
				pedout.print("sample_" + i + " ");
			}
			pedout.print(1 + " ");
			pedout.print(0 + " ");
			pedout.print(0 + " ");
			pedout.print(1 + " ");
			pedout.print(phenotype[i][0] + " ");

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

			if (inbredArgs.isFID())
			{
				phe.print(inbredArgs.getFamIDPrefix() + "sample_" + i + " " + 1);
				breed.println(inbredArgs.getFamIDPrefix() + "sample_" + i + " " + 1 + " " + BV[i]);				
			}
			else
			{
				phe.print("sample_" + i + " " + 1);
				breed.println("sample_" + i + " " + 1 + " " + BV[i]);				
			}
			for (int j = 0; j < rep; j++)
			{
				phe.print(" " + phenotype[i][j]);
			}
			phe.println();
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
		phe.close();
		geno.close();
		breed.close();
	}

	private void writeEffFile()
	{
		PrintWriter eff = null;
		try
		{
			eff = new PrintWriter(new BufferedWriter(new FileWriter(inbredArgs.getOutRoot() + ".rnd")));
		}
		catch (IOException e)
		{
			Logger.handleException(e,
					"An exception occurred when writing files.");
		}

		for (int i = 0; i < M; i++)
		{
			eff.println("rs" + i + " " + A1 + " " + effect[i]);
		}
		eff.close();
	}
	
	private SimulationInbredCommandArguments inbredArgs;

	private RandomDataImpl rnd = new RandomDataImpl();
	private long seed;

	private int M;
	private int nullM;
	private int sample;
	private int rep;

	private double[][] genotype;
	private double[] BV;
	private double[][] phenotype;

	private double[] effect;
	private double[] freq;

	private double[] dprime;
	private double[] LD;

	private double h2;

	private String A1 = "A";
	private String A2 = "C";

}
