package gear.subcommands.metapc;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.EigenDecompositionImpl;
import org.apache.commons.math.stat.StatUtils;

import gear.subcommands.metapc.freader.FReader;
import gear.subcommands.metapc.freader.FStat;
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

		Logger.printUserLog("Principal component analysis for summary statistics.\n");

		Logger.printUserLog(mpcArgs.toString());

		initial();
		calculateGRM();
		EigenAnalysis();
	}

	private void initial()
	{
		boolean[] FileKeep = new boolean[mpcArgs.getMetaFile().length];
		Arrays.fill(FileKeep, true);
		fReader = new FReader(mpcArgs.getMetaFile(), FileKeep,
				mpcArgs.getKeys(), mpcArgs.isQT(), mpcArgs.isGZ(),
				mpcArgs.isChr(), mpcArgs.getChr());

		fReader.Start();

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

		HashMap<String, ArrayList<Float>> snpFrq = NewIt.newHashMap();
		HashMap<String, FStat> SumStat1 = fReader.getMetaStat().get(0);

		HashMap<String, ArrayList<Integer>> snpCntTable = fReader.getMetaSNPTable();

		HashSet<String> atgcSNP = NewIt.newHashSet();
		for (String snp : snpCntTable.keySet())
		{
			if(SumStat1.containsKey(snp))
			{
				FStat ms = SumStat1.get(snp);

				if (SNPMatch.isAmbiguous(ms.getA1(), ms.getA2()))
				{
					atgcSNP.add(snp);
				}					

				ArrayList<Integer> cnt = snpCntTable.get(snp);
				if (cnt.get(cnt.size()-1) == fReader.getNumMetaFile())
				{
					ArrayList<Float> fq = NewIt.newArrayList();
					fq.add(ms.getEffect());
					snpFrq.put(snp, fq);
				}
			}
		}

		for (int i = 1; i < fReader.getNumMetaFile(); i++)
		{
			HashMap<String, FStat> SumStat2 = fReader.getMetaStat().get(i);			
			for (String snp : snpCntTable.keySet())
			{
				ArrayList<Integer> cnt = snpCntTable.get(snp);

				if (cnt.get(cnt.size() -1) == fReader.getNumMetaFile())
				{

					FStat ms1 = SumStat1.get(snp);
					FStat ms2 = SumStat2.get(snp);

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
						atgcSNP.add(snp);
//						snpFrq.remove(snp);
//						continue;
					}

					ArrayList<Float> fq = snpFrq.get(snp);
					if (lineup)
					{
						fq.add(ms2.getEffect());
					}
					else
					{
						if (mpcArgs.isBeta())
						{
							fq.add(-1*ms2.getEffect());
						}
						else
						{
							fq.add(1-ms2.getEffect());							
						}
					}
				}
			}
		}

		if(!mpcArgs.isKeepATGC() )
		{
			if (atgcSNP.size() > 0)
			{
				Logger.printUserLog(atgcSNP.size() + " palindromic marker(s) have been removed.");
				for(Iterator<String> s =  atgcSNP.iterator(); s.hasNext();)
				{
					String snp = s.next();
					snpFrq.remove(snp);
				}
			}
		}

		Logger.printUserLog(snpFrq.size() + " consensus markers have been found and used for analysis.");

		double[][] mg = new double[snpFrq.size()][fReader.getNumMetaFile()];

		int cnt = 0;
		for (String snp:snpFrq.keySet())
		{
			SNPlist.add(snp);
			ArrayList<Float> fq = snpFrq.get(snp);
			for (int j = 0; j < fq.size(); j++)
			{
				mg[cnt][j] = fq.get(j).doubleValue();
			}

			double st = StatUtils.variance(mg[cnt]);
			if (st > 1e-6)
			{
				mg[cnt] = StatUtils.normalize(mg[cnt]);
			}
			else
			{
				for (int j = 0; j < fq.size(); j++)
				{
					mg[cnt][j] = 0;
				}
			}
			cnt++;
		}

		PrintStream frqWriter = FileUtil.CreatePrintStream(new String(mpcArgs.getOutRoot() + ".msnp"));

		for(int i = 0; i < mg.length; i++)
		{
			frqWriter.println(SNPlist.get(i));
		}
		frqWriter.close();

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
		PrintStream grmWriter = FileUtil.CreatePrintStream(new String(mpcArgs.getOutRoot() + ".crm"));
		Array2DRowRealMatrix rm = new Array2DRowRealMatrix(mGRM);

		for (int i = 0; i < rm.getRowDimension(); i++)
		{
			for (int j = 0; j < rm.getColumnDimension(); j++)
			{
				grmWriter.print(rm.getEntry(i, j) + " ");
			}
			grmWriter.println();
		}
		grmWriter.close();

		EigenDecompositionImpl ed = new EigenDecompositionImpl(rm.copy(), 1e-6);

		double[][] ev = new double[rm.getRowDimension()][rm.getColumnDimension()];
		for (int i = 0; i < rm.getColumnDimension(); i++)
		{
			ev[i] = ed.getEigenvector(i).toArray();
		}

		Array2DRowRealMatrix evM = new Array2DRowRealMatrix(ev);
		Array2DRowRealMatrix evMat = (Array2DRowRealMatrix) evM.transpose();


		PrintStream evaWriter = FileUtil.CreatePrintStream(new String(mpcArgs.getOutRoot() + ".mval"));
		double[] eR=ed.getRealEigenvalues();

		for (int i = 0; i < rm.getRowDimension(); i++)
		{
			evaWriter.print(eR[i] + " ");
			evaWriter.println();
		}
		evaWriter.close();

		PrintStream eveWriter = FileUtil.CreatePrintStream(new String(mpcArgs.getOutRoot() + ".mvec"));

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

	private MetaPCCommandArguments mpcArgs;

	private FReader fReader;
	private double[][] mGRM;
	private ArrayList<String> SNPlist = NewIt.newArrayList();
}
