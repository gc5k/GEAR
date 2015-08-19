package gear.subcommands.metapc;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.EigenDecompositionImpl;
import org.apache.commons.math.stat.StatUtils;

import gear.gwassummary.GWASReader;
import gear.gwassummary.MetaStat;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;
import gear.util.SNPMatch;

public class MetaPCCommandImpl extends CommandImpl
{

	@Override
	public void execute(CommandArguments cmdArgs)
	{
		mpcArgs = (MetaPCCommandArguments) cmdArgs;

		if (mpcArgs.isQT())
		{
			Logger.printUserLog("Analysing summary statistics analysis for quantitative traits.\n");
		}
		else
		{
			Logger.printUserLog("Analysing summary statistics analysis for case-contrl studies.\n");
		}

		initial();
		calculateGRM();
		EigenAnalysis();
	}

	private void initial()
	{
		boolean[] FileKeep = new boolean[mpcArgs.getMetaFile().length];
		Arrays.fill(FileKeep, true);
		gReader = new GWASReader(mpcArgs.getMetaFile(), FileKeep,
				mpcArgs.getKeys(), mpcArgs.isQT(), mpcArgs.isGZ(),
				mpcArgs.isChr(), mpcArgs.getChr());

		gReader.Start(mpcArgs.isFrq());

		if (mpcArgs.isFrq())
		{
			Logger.printUserLog("Calculating allele frequency difference, and Fst.");
		}

		int NumMetaFile = mpcArgs.getMetaFile().length;

		if (NumMetaFile < 2)
		{
			Logger.printUserError("At least two summary statistic files should be specified.\n");
			Logger.printUserError("GEAR quitted.\n");
			System.exit(0);
		}

		mGRM = new double[NumMetaFile][NumMetaFile];

		// reading meta files
	}

	private void calculateGRM()
	{

		int cntAmbiguous = 0;
		HashMap<String, ArrayList<Float>> snpFrq = NewIt.newHashMap();
		HashMap<String, MetaStat> SumStat1 = gReader.getMetaStat().get(0);

		HashMap<String, ArrayList<Integer>> snpCntTable = gReader.getMetaSNPTable();

		HashSet<String> badSNP = NewIt.newHashSet();
		for (String snp : snpCntTable.keySet())
		{
			if(SumStat1.containsKey(snp))
			{
				MetaStat ms = SumStat1.get(snp);

				if (SNPMatch.isAmbiguous(ms.getA1(), ms.getA2()))
				{
					badSNP.add(snp);
					cntAmbiguous++;
					continue;
				}

				ArrayList<Integer> cnt = snpCntTable.get(snp);
				if (cnt.get(cnt.size()-1) == gReader.getNumMetaFile())
				{
					ArrayList<Float> fq = NewIt.newArrayList();
					fq.add(ms.getEffect());
					snpFrq.put(snp, fq);
				}
			}
		}

		for (int i = 1; i < gReader.getNumMetaFile(); i++)
		{
			HashMap<String, MetaStat> SumStat2 = gReader.getMetaStat().get(i);			
			for (String snp : snpCntTable.keySet())
			{
				
				if (badSNP.contains(snp))
				{
					continue;
				}
				ArrayList<Integer> cnt = snpCntTable.get(snp);

				if (cnt.get(cnt.size() -1) == gReader.getNumMetaFile())
				{
//					System.out.println(snp);

					MetaStat ms1 = SumStat1.get(snp);
					MetaStat ms2 = SumStat2.get(snp);

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
						badSNP.add(snp);
						snpFrq.remove(snp);
						continue;
					}

					ArrayList<Float> fq = snpFrq.get(snp);
					if (lineup)
					{
						fq.add(ms2.getEffect());
					}
					else
					{
						fq.add(1-ms2.getEffect());
					}
				}
			}
		}

		Logger.printUserLog(cntAmbiguous + " SNPs have been removed.");
		Logger.printUserLog(snpFrq.size() + " consensus markers have been found.");

		double[][] mg = new double[snpFrq.size()][gReader.getNumMetaFile()];

		int cnt = 0;
		for (String snp:snpFrq.keySet())
		{
			SNPlist.add(snp);
			ArrayList<Float> fq = snpFrq.get(snp);
			for (int j = 0; j < fq.size(); j++)
			{
				mg[cnt][j] = fq.get(j);
			}
			mg[cnt] = StatUtils.normalize(mg[cnt]);
			cnt++;
		}

		for (int i = 0; i < mGRM.length; i++)
		{
			for (int j = 0; j < mGRM[i].length; j++)
			{
				for (int k = 0; k < mg.length; k++)
				{
					mGRM[i][j] += mg[k][i] * mg[k][j]; 
				}
				mGRM[i][j] /= mg.length;
			}
		}
	}

	private void EigenAnalysis()
	{
		Array2DRowRealMatrix rm = new Array2DRowRealMatrix(mGRM);
		EigenDecompositionImpl ed = new EigenDecompositionImpl(rm.copy(), 1e-6);

		double[][] ev = new double[rm.getRowDimension()][rm.getColumnDimension()];
		for (int i = 0; i < rm.getColumnDimension(); i++)
		{
			ev[i] = ed.getEigenvector(i).toArray();
		}

		Array2DRowRealMatrix evM = new Array2DRowRealMatrix(ev);
		Array2DRowRealMatrix evMat = (Array2DRowRealMatrix) evM.transpose();

		PrintStream grmWriter = FileUtil.CreatePrintStream(new String(mpcArgs.getOutRoot() + ".mgrm"));
		
		for (int i = 0; i < rm.getRowDimension(); i++)
		{
			for (int j = 0; j < rm.getColumnDimension(); j++)
			{
				grmWriter.print(rm.getEntry(i, j) + " ");
			}
			grmWriter.println();
		}
		grmWriter.close();
		
		PrintStream evaWriter = FileUtil.CreatePrintStream(new String(mpcArgs.getOutRoot() + ".eigenval"));
		double[] eR=ed.getRealEigenvalues();

		for (int i = 0; i < rm.getRowDimension(); i++)
		{
			evaWriter.print(eR[i] + " ");
			evaWriter.println();
		}
		evaWriter.close();

		PrintStream eveWriter = FileUtil.CreatePrintStream(new String(mpcArgs.getOutRoot() + ".eigenvec"));
		
		for (int i = 0; i < evMat.getRowDimension(); i++)
		{
			for (int j = 0; j < evMat.getColumnDimension(); j++)
			{
				eveWriter.print(evMat.getEntry(i, j) + " ");
			}
			eveWriter.println();
		}
		eveWriter.close();
	}
	
	Array2DRowRealMatrix rm;
	private MetaPCCommandArguments mpcArgs;

	private GWASReader gReader;
	private double[][] mGRM;
	private ArrayList<String> SNPlist = NewIt.newArrayList();

}
