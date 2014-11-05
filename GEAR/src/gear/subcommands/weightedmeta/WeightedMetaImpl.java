package gear.subcommands.weightedmeta;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.Set;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistributionImpl;

import gear.gwassummary.GWASReader;
import gear.gwassummary.MetaStat;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.BufferedReader;
import gear.util.Logger;
import gear.util.NewIt;
import gear.util.SNPMatch;

public class WeightedMetaImpl extends CommandImpl
{

	@Override
	public void execute(CommandArguments cmdArgs)
	{
		wMetaArgs = (WeightedMetaArguments) cmdArgs;

		if (wMetaArgs.isQT())
		{
			Logger.printUserLog("Analysing summary statistics analysis for quantitative traits.\n");			
		}
		else
		{
			Logger.printUserLog("Analysing summary statistics analysis for case-contrl studies.\n");			
		}

		gReader = new GWASReader(wMetaArgs.getMetaFile(), wMetaArgs.getKeys(), wMetaArgs.isQT(), wMetaArgs.isGZ());
		
		if (gReader.getNumMetaFile() < 2)
		{
			Logger.printUserError("At least two summary statistic files should be specified.\n");
			Logger.printUserError("GEAR quitted.\n");
			System.exit(0);
		}

		generateCorMatrix();			
		MetaAnalysis();
	}

	private void generateCorMatrix()
	{
		corMat = new double[gReader.getNumMetaFile()][gReader.getNumMetaFile()];
		if (wMetaArgs.getCMFile() == null)
		{
			for(int i = 0; i < corMat.length; i++)
			{
				corMat[i][i] = 1;
			}
			Logger.printUserLog("No correlation matrix is specified. The default correlation (digonal matrix) will be used.");
		}
		else
		{
			int NumMetaFile = gReader.getNumMetaFile();
			BufferedReader bf = BufferedReader.openTextFile(wMetaArgs.getCMFile(), "cm file.");
			Logger.printUserLog("Reading '" + wMetaArgs.getCMFile() + "'.");

			String[] d = null;
			int cnt = 0;
			while ( (d = bf.readTokens())!= null )
			{
				if (d.length != NumMetaFile )
				{
					Logger.printUserError("incorrect '" + wMetaArgs.getCMFile() + "'.");
					System.exit(0);
				}
				for(int i = 0; i < d.length; i++)
				{
					corMat[cnt][i] = Double.parseDouble(d[i]);
				}
				cnt++;
			}
			for(int i = 0; i < corMat.length; i++)
			{
				for(int j = 0; j < i; j++)
				{
					corMat[j][i] = corMat[i][j];
				}
			}
			Logger.printUserLog(corMat.length + "X" + corMat.length + " correlation matrix has been read in.");
		}
	}

	private void MetaAnalysis()
	{
		Logger.printUserLog("Starting meta-analysis...");
		int cnt = 0;
		int atgcCnt = 0;
		Set<String> keys = gReader.getMetaSNPTable().keySet();
		for (Iterator<String> e=keys.iterator(); e.hasNext();)
		{
			cnt++;
			String key = e.next();
			System.out.println(key);
			ArrayList<Integer> Int = gReader.getMetaSNPTable().get(key);
			CovMatrix covMat = new CovMatrix(key, Int, corMat, gReader, wMetaArgs.getGC(), wMetaArgs.getGCInflationOnly());
			GMRes gr = MetaSNP(covMat);
			if (gr.getIsAmbiguous())
			{
				atgcCnt++;
				if (!wMetaArgs.isKeepATGC())
				{
					continue;
				}
			}
			grArray.add(gr);
		}
		Collections.sort(grArray);
		Logger.printUserLog("In total "+ cnt + " loci have been analyzed for meta-analysis.");
		if (!wMetaArgs.isKeepATGC())
		{
			Logger.printUserLog(atgcCnt + " ambiguous loci have been eliminated.");
		}
		PrintGMresults();
	}

	private GMRes MetaSNP(CovMatrix covMat)
	{
		String SNP = covMat.getSNP();
		int[] idx = covMat.getCohortIdx();
		int cohort = idx.length;
		double[] Weight = covMat.getWeights();
		double gse = covMat.getGSE();

		StringBuffer direct = new StringBuffer();

		for(int i = 0; i < gReader.getNumMetaFile(); i++)
		{
			direct.append('?');
		}
		GMRes gr = new GMRes(cohort);

		double gb = 0;

		MetaStat ms = null;
		boolean isAmbiguousLocus = false;
		for (int i = 0; i < idx.length; i++)
		{
			char sign = '+';
			double b = 0;
			boolean match = true;
			ms = gReader.getMetaStat().get(idx[i]).get(SNP);
			if (i == 0)
			{
				b = ms.getEffect();
				gr.SetSNP(ms.getSNP());
				gr.SetChr(ms.getChr());
				gr.SetBP(ms.getBP());
				gr.SetA1(ms.getA1());
				gr.SetA2(ms.getA2());
				isAmbiguousLocus = SNPMatch.isAmbiguous(ms.getA1(), ms.getA2());
			}
			else
			{
				match = SNPMatch.isAllelsMatchForTwoLoci(gr.GetA1(), gr.GetA2(), ms.getA1(), ms.getA2());
				if (match)
				{
					if (gr.GetA1() == ms.getA1() || gr.GetA1() == SNPMatch.Flip(ms.getA1())) //match A1 in the second meta
					{
						b = ms.getEffect();
					}
					else if (gr.GetA1() == ms.getA2() || gr.GetA1() == SNPMatch.Flip(ms.getA2())) //match A2 in the second meta
					{
						b = -1 * ms.getEffect();
					}
				}
				else
				{
					sign = ',';
				}
			}

			if(match)
			{
				gb += b * Weight[i];
				sign = b > 0 ? '+':'-';
			}
			direct.setCharAt(idx[i], sign);
		}
		double z = gb/gse;
		double p = 1;
		try
		{
			p = (1-unitNormal.cumulativeProbability(Math.abs(z)))*2;
		}
		catch (MathException e)
		{
			Logger.printUserError(e.toString());
		}
		gr.SetAmbi(isAmbiguousLocus);
		gr.SetB(gb);
		gr.SetSE(gse);
		gr.SetZ(z);
		gr.SetP(p);
		gr.SetDirect(direct.toString());
		return gr;
	}

	private void PrintGMresults()
	{
        PrintWriter writer = null;
        try
        {
        	writer = new PrintWriter(new BufferedWriter(new FileWriter(wMetaArgs.getOutRoot()+".gmeta")));
        	Logger.printUserLog("Writting detailed test statistics into '"+wMetaArgs.getOutRoot() + ".gmeta.'\n");
        }
		catch (IOException e)
		{
			Logger.handleException(e, "An I/O exception occurred when writing '" + wMetaArgs.getOutRoot() + ".gmeta" + "'.\n");
		}

		for(int i = 0; i < grArray.size(); i++)
		{
			GMRes gr = grArray.get(i);
			if(i == 0)
			{
				writer.write(gr.printTitle() + "\n");
			}
			writer.write(gr.toString()+ "\n");
		}
		writer.close();
	}

	private WeightedMetaArguments wMetaArgs;
	private GWASReader gReader;
	private double[][] corMat;
	private NormalDistributionImpl unitNormal = new NormalDistributionImpl(0, 1);
	private ArrayList<GMRes> grArray = NewIt.newArrayList();
}
