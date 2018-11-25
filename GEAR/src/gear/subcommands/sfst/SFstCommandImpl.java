package gear.subcommands.sfst;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;

import gear.gwassummary.GWASReader;
import gear.gwassummary.MetaStat;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;
import gear.util.SNPMatch;

public class SFstCommandImpl extends CommandImpl
{

	@Override
	public void execute(CommandArguments cmdArgs)
	{
		sfstArgs = (SFstCommandArguments) cmdArgs;

		if (sfstArgs.isQT())
		{
			Logger.printUserLog("Analysing summary statistics analysis for quantitative traits.\n");
		}
		else
		{
			Logger.printUserLog("Analysing summary statistics analysis for case-contrl studies.\n");
		}

		initial();

		// generating matrix
		String[] MetaFile = gReader.getMetaFile();
		int FileSize = sfstArgs.getTop() > 0 ? sfstArgs.getTop()
				: (MetaFile.length - 1);
		for (int i = 0; i < FileSize; i++)
		{
			for (int j = (i + 1); j < MetaFile.length; j++)
			{
				Logger.printUserLog("");
				Logger.printUserLog("File pair: " + (i + 1) + "-" + (j + 1));
				if (sfstArgs.isQT())
				{
					double[] size = sfstArgs.getQTsize();
					Kappa = 2 / (Math.sqrt(size[i] / size[j]) + Math
							.sqrt(size[j] / size[i]));
					Logger.printUserLog("Sample sizes for '" + MetaFile[i] + "': " + size[i]);
					Logger.printUserLog("Sample sizes for '" + MetaFile[j] + "': " + size[j]);
				}
				else
				{
					double[] size = sfstArgs.getCCsize();
					R1 = size[i * 2] / size[i * 2 + 1];
					R2 = size[j * 2] / size[j * 2 + 1];
					double s1 = size[i * 2] + size[i * 2 + 1];
					double s2 = size[j * 2] + size[j * 2 + 1];
					Kappa = 2 / (Math.sqrt(s1 / s2) + Math.sqrt(s2 / s1));
					Logger.printUserLog("Sample size for '" + MetaFile[i] + "': " + size[i * 2] + " cases, " + size[i * 2 + 1] + " controls; R1 = " + R1 + ".");
					Logger.printUserLog("Sample size for '" + MetaFile[j] + "': " + size[j * 2] + " cases, " + size[j * 2 + 1] + " controls; R2 = " + R2 + ".");
				}
				Logger.printUserLog("Kappa: " + Kappa);

				calculateFst(i, j);
			}
		}
		WriteFstMat();
		Logger.printUserLog("=========================================================");
		Logger.printUserLog("Results have been saved in '" + sfstArgs
				.getOutRoot() + ".lmat'.");
	}

	private void initial()
	{
		boolean[] FileKeep = new boolean[sfstArgs.getMetaFile().length];
		Arrays.fill(FileKeep, true);
		gReader = new GWASReader(sfstArgs.getMetaFile(), FileKeep,
				sfstArgs.getKeys(), sfstArgs.isQT(), sfstArgs.isGZ(),
				sfstArgs.getChr() != null? true:false,
				sfstArgs.getChr() != null?sfstArgs.getChr():sfstArgs.getNotChr());

		gReader.Start(sfstArgs.isFrq());

		Me = sfstArgs.getMe();
		if (Me > 0)
		{
			Logger.printUserLog("Set the effective number of marker: " + Me);
		}
		else
		{
			Logger.printUserLog("Using all markers.");
		}

		if (sfstArgs.isFrq())
		{
			Logger.printUserLog("Calculating allele frequency difference, and Fst.");
		}

		int NumMetaFile = sfstArgs.getMetaFile().length;

		if (NumMetaFile < 2)
		{
			Logger.printUserLog("At least two summary statistic files should be specified. GEAR quitted.");
			System.exit(0);
		}

		fstMat = new double[NumMetaFile][NumMetaFile];

		// reading meta files
	}

	private void calculateFst(int idx1, int idx2)
	{
		ArrayList<FstUnit> FstArray = NewIt.newArrayList();

		// DescriptiveStatistics T0 = new DescriptiveStatistics();

		int cntAmbiguous = 0;
		HashMap<String, MetaStat> SumStat1 = gReader.getMetaStat().get(idx1);
		HashMap<String, MetaStat> SumStat2 = gReader.getMetaStat().get(idx2);

		ArrayList<String> snpArray = gReader.getMetaSNPArray().get(idx1);

		int[][] KeyIdx = gReader.getKeyIndex();
		for (String snp : snpArray)
		{
			if (!SumStat2.containsKey(snp))
			{
				continue;
			}
			MetaStat ms1 = SumStat1.get(snp);
			MetaStat ms2 = SumStat2.get(snp);

			if (KeyIdx[idx1][GWASReader.SE] != -1)
			{
				if (SNPMatch.isAmbiguous(ms1.getA1(), ms1.getA2()))
				{
					cntAmbiguous++;
					continue;
				}
			}
			if (KeyIdx[idx2][GWASReader.SE] != -1)
			{
				if (SNPMatch.isAmbiguous(ms2.getA1(), ms2.getA2()))
				{
					cntAmbiguous++;
					continue;
				}
			}

			boolean lineup = true;
			if (ms1.getA1() == ms2.getA1() || ms1.getA1() == SNPMatch.Flip(ms2
					.getA1())) // match A1 in the second meta
			{
			}
			else if (ms1.getA1() == ms2.getA2() || ms1.getA1() == SNPMatch
					.Flip(ms2.getA2())) // match A2 in the second meta
			{
				lineup = false;
			}
			else
			{
				cntAmbiguous++;
				continue;
			}

			double s1, s2;
			if (sfstArgs.isNoWeight())
			{
				s1 = sfstArgs.getNe();
				s2 = sfstArgs.getNe();
			}
			else
			{
				if (sfstArgs.isQT())
				{
					s1 = sfstArgs.getQTsize()[idx1];
					s2 = sfstArgs.getQTsize()[idx2];
				}
				else
				{
					s1 = sfstArgs.getCCsize()[idx1 * 2] + sfstArgs.getCCsize()[idx1 * 2 + 1];
					s2 = sfstArgs.getCCsize()[idx2 * 2] + sfstArgs.getCCsize()[idx2 * 2 + 1];
				}
			}

			FstArray.add(new FstUnit(ms1, ms2, lineup, s1, s2));
		}

		if (cntAmbiguous > 0)
		{
			if (cntAmbiguous == 1)
			{
				Logger.printUserLog("Removed " + cntAmbiguous + " ambiguous locus (AT/GC).");
			}
			else
			{
				Logger.printUserLog("Removed " + cntAmbiguous + " ambiguous loci (AT/GC).");
			}
		}
		Logger.printUserLog("Found " + FstArray.size() + " consensus summary statistics between these two files.");

		if (FstArray.size() < ((int) (Me * sfstArgs.getMeFrac())))
		{
			Logger.printUserLog("Too few overlapping snps, skip this pair of files.");
			return;
		}

		// select independent snps
		Collections.sort(FstArray);

		int[] selIdx = null;
		if (Me < 0)
		{// use all
			selIdx = new int[FstArray.size()];
			for (int i = 0; i < selIdx.length; i++)
				selIdx[i] = i;
		}
		else if (FstArray.size() <= Me)
		{// use available ones
			selIdx = new int[FstArray.size()];
			for (int i = 0; i < selIdx.length; i++)
				selIdx[i] = i;
		}
		else
		{// use Me
			selIdx = new int[(int) Math.ceil(Me)];
			for (int i = 0; i < Me; i++)
				selIdx[i] = (int) Math.floor((i * 1.0 + 1) / Me * FstArray
						.size()) - 1;
		}

		double fst = 0;

		for (int i = 0; i < selIdx.length; i++)
		{
			FstUnit lu = FstArray.get(selIdx[i]);
			fst += lu.getFstBW()/selIdx.length;
		}
		FstDistance(FstArray, selIdx);

		Logger.printUserLog("Fst is " + fst);
//		Logger.printUserLog("FstD is " + FstDistance);
		fstMat[idx2][idx1] = fstMat[idx1][idx2] = fst;

		if (sfstArgs.isVerbose())
		{
			VerboseGZ(FstArray, idx1, idx2, selIdx);
		}
	}

	private void VerboseGZ(ArrayList<FstUnit> FstArray, int idx1, int idx2, int[] selIdx)
	{
		BufferedWriter GZ = FileUtil
				.ZipFileWriter(new String(
						sfstArgs.getOutRoot() + "." + (idx1 + 1) + "-" + (idx2 + 1) + ".sFst.gz"));

		try
		{
			GZ.append(titleLine);
		}
		catch (IOException e)
		{
			Logger.handleException(e,
									"error in writing " + new String(
											sfstArgs.getOutRoot() + "." + (idx1 + 1) + "-" + (idx2 + 1) + ".sFst.gz"));
		}

		for (int i = 0; i < selIdx.length; i++)
		{
			FstUnit lu = FstArray.get(selIdx[i]);
			MetaStat ms1 = lu.getMetaStat1();
			MetaStat ms2 = lu.getMetaStat2();
			try
			{
				GZ.write(ms1.getSNP() + "\t" + ms1.getChr() + "\t" + ms1
						.getBP() + "\t" + ms1.getA1() + "\t" + lu.getB1() + "\t" + ms1
						.getSE() + "\t" + ms1.getP() + "\t" + lu.getB2() + "\t" + ms2
						.getSE() + "\t" + ms2.getP() + "\t" + lu.getFstBW() + "\n");
			}
			catch (IOException e)
			{
				Logger.handleException(e,
										"error in writing " + new String(
												sfstArgs.getOutRoot() + "." + (idx1 + 1) + "-" + (idx2 + 1) + ".sFst.gz"));
			}	
		}
		
		try
		{
			GZ.close();
		}
		catch (IOException e)
		{
			Logger.handleException(e,
									"error in writing " + new String(
											sfstArgs.getOutRoot() + "." + (idx1 + 1) + "-" + (idx2 + 1) + ".sFst.gz"));
		}

		Logger.printUserLog("Save the fst results between this pair of files into " + sfstArgs.getOutRoot() + "." + (idx1 + 1) + "-" + (idx2 + 1) + ".sFst.gz");
	}

	private void WriteFstMat()
	{
		PrintStream FstWriter = FileUtil.CreatePrintStream(new String(sfstArgs.getOutRoot() + ".fst"));
		for (int i = 0; i < fstMat.length; i++)
		{
			for (int j = 0; j < fstMat[i].length; j++)
			{
				FstWriter.print(String.format("%.5f", fstMat[i][j]) + " ");
			}
			FstWriter.println();
		}
		FstWriter.close();
	}

	private void FstDistance(ArrayList<FstUnit> fstUnit, int[] idx)
	{
		double fstD1 = 0;
		double fstD2 = 0;
//		double fstD3 = 0;
		double beta1_n = 0;
		double beta1_d = 0;
		double beta2_n = 0;
		double beta2_d = 0;

		for(int i = 0; i < idx.length; i++)
		{
			FstUnit fu = fstUnit.get(idx[i]);
			double p1 = fu.getB1();
			double p2 = fu.getB2();
			double n1 = fu.getN1() * 2;
			double n2 = fu.getN2() * 2;
			
			beta1_n += (n1*2)/(n1*2-1)*(1-p1*p1 - (1-p1) * (1-p1));
			beta1_d += (1-p1*p2 - (1-p1) * (1-p2));
			
			beta2_n += (n2*2)/(n2*2-1)*(1-p2*p2 - (1-p2) * (1-p2));
			beta2_d += (1-p1*p2 - (1-p1) * (1-p2));
			
		}
		fstD1 = 1 - beta1_n / beta1_d;
		fstD2 = 1 - beta2_n / beta2_d;

		FstDistance = (fstD1 + fstD2)/2;
	}

	private double Me = 30000;
	private double R1 = 1;
	private double R2 = 1;
	private double Kappa = 1;
	private SFstCommandArguments sfstArgs;

	private String titleLine = "SNP\tChr\tBp\tA1\tFrq1\tSE1\tP1\tFrq2\tSE2\tP2\tFst\n";

	private GWASReader gReader;
	private double FstDistance = 0;

	private double[][] fstMat;

}
