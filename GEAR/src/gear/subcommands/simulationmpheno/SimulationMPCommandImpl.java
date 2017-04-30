package gear.subcommands.simulationmpheno;

import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.random.RandomDataImpl;
import org.apache.commons.math.stat.StatUtils;
import org.apache.commons.math.linear.CholeskyDecompositionImpl;
import org.apache.commons.math.linear.NonSquareMatrixException;
import org.apache.commons.math.linear.NotPositiveDefiniteMatrixException;
import org.apache.commons.math.linear.NotSymmetricMatrixException;

import gear.ConstValues;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.BufferedReader;
import gear.util.Logger;
import gear.util.Sample;
import gear.util.pop.PopStat;

public class SimulationMPCommandImpl extends CommandImpl
{

	@Override
	public void execute(CommandArguments cmdArgs)
	{
		mpArgs = (SimulationMPCommandArguments)cmdArgs;

		sample = mpArgs.getSampleSize();
		h2 = mpArgs.getHsq();
		M = mpArgs.getMarkerNum();
		nullM = mpArgs.getNullMarkerNum();
		seed = mpArgs.getSeed();
		rnd.reSeed(seed);
		rep = mpArgs.getRep();

		readCorMatrix();
		getFreq();
		getEffect();
		getDPrime();
		calLD();
		generateSampleNoSelection();
		
		if (mpArgs.isMakeBed())
		{
			writeBFile();
		} 
		else
		{
			writeFile();
		}
	}

	private void readCorMatrix()
	{
		corMat = new double[h2.length][h2.length];
		corEMat = new double[h2.length][h2.length];
		for(int i = 0; i < corMat.length; i++)
		{
			corMat[i][i] = 1;
			corEMat[i][i] = 1;
		}

		if (mpArgs.getCMFile() != null)
		{
			BufferedReader bf = BufferedReader.openTextFile(mpArgs.getCMFile(), "cm file.");
			Logger.printUserLog("Reading '" + mpArgs.getCMFile() + "'.");

			String[] d = null;
			int cIdx = 0;
			while ( (d = bf.readTokens())!= null )
			{
				if (d.length != h2.length)
				{
					Logger.printUserError("incorrect '" + mpArgs.getCMFile() + "'.");
					System.exit(0);
				}
				int c = 0;
				for(int i = 0; i < d.length; i++)
				{
					corMat[cIdx][c++] = Double.parseDouble(d[i]);
				}
				cIdx++;
			}
			Logger.printUserLog(corMat.length + "X" + corMat.length + " correlation matrix has been read in.");	
		}

		if (mpArgs.getCMEFile() != null)
		{
			BufferedReader bf = BufferedReader.openTextFile(mpArgs.getCMEFile(), "cm file.");
			Logger.printUserLog("Reading '" + mpArgs.getCMEFile() + "'.");

			String[] d = null;
			int cIdx = 0;
			while ( (d = bf.readTokens())!= null )
			{
				if (d.length != h2.length)
				{
					Logger.printUserError("incorrect '" + mpArgs.getCMEFile() + "'.");
					System.exit(0);
				}
				int c = 0;
				for(int i = 0; i < d.length; i++)
				{
					corEMat[cIdx][c++] = Double.parseDouble(d[i]);
				}
				cIdx++;
			}
			Logger.printUserLog(corMat.length + "X" + corMat.length + " environmental correlation matrix has been read in.");	
		}

		CholeskyDecompositionImpl cholCImpl = null;
		CholeskyDecompositionImpl cholCEImpl = null;

		try {
			cholCImpl = new CholeskyDecompositionImpl(new Array2DRowRealMatrix(corMat));
			cholCEImpl = new CholeskyDecompositionImpl(new Array2DRowRealMatrix(corEMat));
		} catch (NonSquareMatrixException e) {
			e.printStackTrace();
		} catch (NotSymmetricMatrixException e) {
			e.printStackTrace();
		} catch (NotPositiveDefiniteMatrixException e) {
			e.printStackTrace();
		} catch (IllegalArgumentException e) {
			e.printStackTrace();
		} catch (NullPointerException e) {
			e.printStackTrace();
		}
		RealMatrix L = cholCImpl.getL();
		System.out.println(L);
		cholC = L.getData();

		RealMatrix L1 = cholCEImpl.getL();
		System.out.println(L1);
		cholCE = L1.getData();

	}

	private void generateSampleNoSelection()
	{
		phenotype = new double[sample][h2.length][rep];
		BV = new double[sample][h2.length];
		genotype = new double[sample][M];

		for (int i = 0; i < sample; i++)
		{
			RealMatrix chr = SampleChromosome();
			genotype[i] = chr.getColumn(0);

			for (int t = 0; t < h2.length; t++)
			{
				if (h2[t] == 0)
				{
					for(int j = 0; j < sample; j++)
					{
						BV[i][t] = 0;
					}
				}
				else
				{
					double[] eff = new double[M];
					for (int j = 0; j < M; j++)
					{
						eff[j] = effect[j][t];
					}
					RealMatrix Meffect = new Array2DRowRealMatrix(eff);

					RealMatrix genoEff = chr.transpose().multiply(Meffect);
					double bv = genoEff.getEntry(0, 0);
					BV[i][t] = bv;					
				}
			}
		}

		for (int r = 0; r < rep; r++)
		{
			double[][] res = new double[sample][h2.length];
			for (int j1 = 0; j1 < sample; j1++)
			{
				double[] ze = new double[h2.length];
				for (int j2 = 0; j2 < h2.length; j2++)
				{
					ze[j2] = rnd.nextGaussian(0, 1);
				}

				for (int j2 = 0; j2 < h2.length; j2++)
				{
					for (int j3 = 0; j3 < h2.length; j3++)
					{
						res[j1][j2] += ze[j3] * cholCE[j2][j3];						
					}
				}
			}

			for (int t = 0; t < h2.length; t++)
			{
				double vg = StatUtils.variance((new Array2DRowRealMatrix(BV)).getColumn(t));
				//rescale the phenotype to get the heritability and residual
				double ve = h2[t] == 0 ? 1:vg * (1 - h2[t]) / h2[t];
				double E = Math.sqrt(ve);
				Logger.printUserLog("Trait "+(t+1) + " Vg=" + vg);
				double[] pv=new double[sample];

				for (int i = 0; i < sample; i++)
				{
					phenotype[i][t][r] = BV[i][t] + res[i][t] * E;
					pv[i] = phenotype[i][t][r];
				}
				double Vp=StatUtils.variance(pv);
				Logger.printUserLog("Vp=" + Vp + "; hsq=" + vg/Vp + " for replicate " + (r+1));

			}
		}

//		
//		for (int t = 0; t < h2.length; t++)
//		{
//
//			double vg = StatUtils.variance((new Array2DRowRealMatrix(BV)).getColumn(t));
//			//rescale the phenotype to get the heritability and residual
//			double ve = h2[t] == 0 ? 1:vg * (1 - h2[t]) / h2[t];
//			double E = Math.sqrt(ve);
//			Logger.printUserLog("Vg=" + vg);
//			for (int i = 0; i < rep; i++)
//			{
//				double[] pv=new double[sample];
//				for (int j = 0; j < sample; j++)
//				{
//					phenotype[j][t][i] = BV[j][t] + rnd.nextGaussian(0, E);
//					pv[j] = phenotype[j][t][i];
//				}
//				double Vp=StatUtils.variance(pv);
//				Logger.printUserLog("Vp=" + Vp + "; hsq=" + vg/Vp + " for replicate " + (i+1));
//			}
//		}
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

		if (mpArgs.isPlainFreq())
		{
			Arrays.fill(freq, mpArgs.getFreq());
		}
		else if (mpArgs.isUnifFreq())
		{
			for(int i = 0; i < M; i++)
			{
				freq[i] = rnd.nextUniform(mpArgs.getFreqRangeLow(), mpArgs.getFreqRangeHigh());
			}
		}
		else if (mpArgs.isFreqFile())
		{
			BufferedReader reader = BufferedReader.openTextFile(mpArgs.getFreqFile(), "freq file");
			int c = 0;
			String line = null;
			while ((line = reader.readLine()) != null)
			{
				if(c >= M)
				{
					Logger.printUserLog("Have already read " + M + " allelic frequencies.  Ignore the rest of the content in '" + mpArgs.getFreqFile() + "'.");
				}
				line.trim();
				String[] l = line.split(ConstValues.WHITESPACE_DELIMITER);
				if (l.length < 1) continue;
				freq[c++] = Double.parseDouble(l[0]);
			}
			reader.close();			
		}
	}

	private void getEffect()
	{
		effect = new double[M][h2.length];
		Sample.setSeed(seed);

		int[] idx = Sample.SampleIndex(0, M-1, M-nullM);
		Arrays.sort(idx);

		if (mpArgs.isPlainEffect())
		{
			for (int i = 0; i < h2.length; i++)
			{
				for (int j = 0; j < idx.length; j++) effect[idx[j]][i] = mpArgs.getPolyEffect();				
			}
		}
		else if (mpArgs.isPolyEffect())
		{
			double[][] z = new double[M][h2.length];
			for (int i = 0; i < h2.length; i++)
			{
				for (int j = 0; j < idx.length; j++)
				{
					z[idx[j]][i] = rnd.nextGaussian(0, Math.sqrt(1/( idx.length * 1.0)));
				}
			}

			for (int i = 0; i < h2.length; i++)
			{
				for (int j = 0; j < idx.length; j++)
				{
					for (int k = 0; k < h2.length; k++)
					{
						effect[idx[j]][i] += z[idx[j]][k] * cholC[i][k];						
					}
				}
			}
		}
		else if (mpArgs.isPolyEffectFile())
		{
			BufferedReader reader = BufferedReader.openTextFile(mpArgs.getPolyEffectFile(), "poly effect file.");

			int c = 0;
			String line = null;
			while ((line = reader.readLine()) != null)
			{
				if (c >= M)
				{
					Logger.printUserLog("Have already read " + M + " allelic effects.  Ignore the rest of the content in '" + mpArgs.getFreqFile() + "'.");
					break;
				}

				line.trim();
				String[] l = line.split(ConstValues.WHITESPACE_DELIMITER);
				if (l.length < h2.length) continue;
				for (int i = 0; i < h2.length; i++)
				{
					effect[c][i] = Double.parseDouble(l[i]);						
				}
				c++;
			}
			reader.close();
		}

		for(int i = 0; i < h2.length; i++)
		{
			if (h2[i] == 0)
			{
				for(int j = 0; j < M; j++)
				{
					effect[j][i] = 0;
				}
			}			
		}
	}

	private void getDPrime()
	{
		dprime = new double[M-1];
		if(mpArgs.isPlainLD())
		{
			Arrays.fill(dprime, mpArgs.getLD());
		}
		else if (mpArgs.isRandLD())
		{
			for(int i = 0; i < dprime.length; i++)
			{
				dprime[i] = rnd.nextUniform(mpArgs.getFreqRangeLow(), mpArgs.getFreqRangeHigh());
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
		PrintWriter geno = null;

		try
		{
			bedout = new DataOutputStream(new FileOutputStream(mpArgs.getOutRoot() + ".bed"));

			fam = new PrintWriter(new BufferedWriter(new FileWriter(mpArgs.getOutRoot()
					+ ".fam")));
			bim = new PrintWriter(new BufferedWriter(new FileWriter(mpArgs.getOutRoot()
					+ ".bim")));

			geno = new PrintWriter(new BufferedWriter(new FileWriter(mpArgs.getOutRoot()
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

			fam.println(phenotype[i][0][0] + " ");

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
		bim.close();
		fam.close();


		for (int t = 0; t < h2.length; t++)
		{
			PrintWriter phe = null;
			PrintWriter eff = null;
			PrintWriter breed = null;

			try
			{
				phe = new PrintWriter(new BufferedWriter(new FileWriter(mpArgs.getOutRoot()
						+ "." + (t+1)+ ".phe")));
				eff = new PrintWriter(new BufferedWriter(new FileWriter(mpArgs.getOutRoot() + "." + (t+1) + ".rnd")));
				breed = new PrintWriter(new BufferedWriter(new FileWriter(mpArgs.getOutRoot() + "." + (t+1) +".breed")));
			}
			catch (IOException e)
			{
				e.printStackTrace();
			}

			for (int i = 0; i < sample; i++)
			{
				phe.print("sample_" + i + " " + 1);
				breed.println("sample_" + i + " " + 1 + " " + BV[i][t]);

				for (int j = 0; j < rep; j++)
				{
					phe.print(" " + phenotype[i][t][j]);
				}
				phe.println();
			}

			for (int i = 0; i < M; i++)
			{
				eff.println("rs" + i + " " + A1 + " " + effect[i][t]);
			}

			phe.close();
			eff.close();
			breed.close();
		}

	}

	public void writeFile()
	{
		PrintWriter geno = null;
		PrintWriter map = null;
		PrintWriter pedout = null;

		try
		{
			geno = new PrintWriter(new BufferedWriter(new FileWriter(mpArgs.getOutRoot()
					+ ".add")));
			pedout = new PrintWriter(new BufferedWriter(new FileWriter(mpArgs.getOutRoot()
					+ ".ped")));
			map = new PrintWriter(new BufferedWriter(new FileWriter(mpArgs.getOutRoot()
					+ ".map")));

		}
		catch (IOException e)
		{
			Logger.handleException(e,
					"An exception occurred when writing files.");
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
		geno.close();

		for(int t = 0; t < h2.length; t++)
		{
			PrintWriter phe = null;
			PrintWriter eff = null;
			PrintWriter breed = null;
			try
			{
				phe = new PrintWriter(new BufferedWriter(new FileWriter(mpArgs.getOutRoot()
						+ "." + (t+1) +".phe")));
				eff = new PrintWriter(new BufferedWriter(new FileWriter(mpArgs.getOutRoot() + "." + (t+1) + ".rnd")));
				breed = new PrintWriter(new BufferedWriter(new FileWriter(mpArgs.getOutRoot() + "." + (t+1) + ".breed")));
			}
			catch (IOException e)
			{
				Logger.handleException(e,
						"An exception occurred when writing files.");
			}

			for (int i = 0; i < sample; i++)
			{
				breed.println("sample_" + i + " " + 1 + " " + BV[i][t]);
				phe.print("sample_" + i + " " + 1);
				for (int j = 0; j < rep; j++)
				{
					phe.print(" " + phenotype[i][t][j]);
				}
				phe.println();
			}

			for (int i = 0; i < M; i++)
			{
				eff.println("rs" + i + " " + A1 + " " + effect[i][t]);
			}

			phe.close();
			eff.close();
			breed.close();
		}
	}
	

	private SimulationMPCommandArguments mpArgs;
	
	private RandomDataImpl rnd = new RandomDataImpl();
	private long seed;

	private int M;
	private int nullM;
	private int sample;
	private int rep;

	private double[][] corEMat;
	private double[][] cholCE;
	private double[][] corMat;
	private double[][] cholC;
	private double[][] genotype;
	private double[][] BV;
	private double[][][] phenotype;

	private double[][] effect;
	private double[] freq;

	private double[] dprime;
	private double[] LD;

	private double[] h2;

	private String A1 = "A";
	private String A2 = "C";

}
