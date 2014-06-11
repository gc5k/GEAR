package gear.subcommands.lambdaD;

import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;

import gear.ConstValues;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.BufferedReader;
import gear.util.Logger;
import gear.util.NewIt;
import gear.util.SNPMatch;

public class LambdaDCommandImpl extends CommandImpl
{

	@Override
	public void execute(CommandArguments cmdArgs)
	{
		lamArgs = (LambdaDCommandArguments) cmdArgs;
		initial();
		SumStat1 = readMeta(0);
		SumStat2 = readMeta(1);

		calculateLambdaD();
		Logger.printUserLog("LambdaD median: " + LambdaMedian);
		Logger.printUserLog("Estimated overlapping samples (lambdaD median): " + OSMedian);
		Logger.printUserLog("Estimated rho (lambdaD median): " + rhoMedian);
		if(!lamArgs.isQT())
		{
			Logger.printUserLog("Estiamted overlapping controls: " + OSCtrlMedian);
		}
		Logger.printUserLog("LambdaD mean: " + LambdaMean);
		Logger.printUserLog("Estimated overlapping samples (lambdaD mean): " + OSMean);
		Logger.printUserLog("Estimated rho (lambdaD mean): " + rhoMean);
		if(!lamArgs.isQT())
		{
			Logger.printUserLog("Estiamted overlapping controls: " + OSCtrlMean);			
		}

		printOut();
	}
	
	private void printOut()
	{
		PrintWriter writer = null;
		try 
		{
			writer = new PrintWriter(new BufferedWriter(new FileWriter(lamArgs.getOutRoot() + ".lam")));
		}
		catch (IOException e)
		{
			Logger.handleException(e, "An I/O exception occurred when writing '" + lamArgs.getOutRoot() + ".lam" + "'.");
		}
		writer.println("Kappa: " + Kappa);
		writer.println("LambdaD (median): " + LambdaMedian);
		writer.println("rho (based on LambdaD median: " + rhoMedian);
		writer.println("Overlapping samples: " + OSMedian);
		if(!lamArgs.isQT())
		{
			writer.println("Overlapping controls: " + OSCtrlMedian);
		}
		writer.println();
		writer.println("LambdaD (mean): " + LambdaMean);
		writer.println("rho (based on LambdaD mean): " + rhoMean);
		writer.println("Overlapping samples: " + OSMean);
		if (!lamArgs.isQT())
		{
			writer.println("Overlapping controls: " + OSMedian);
		}
		writer.close();
	}

	private void initial()
	{
		MetaFile[0] = lamArgs.getMeta1();
		MetaFile[1] = lamArgs.getMeta2();

		if (lamArgs.isQT())
		{
			Logger.printUserLog("Summary statistics analysis for quantitative traits.");
			double[] size = lamArgs.getQTsize();
			Kappa = 2 / ( Math.sqrt(size[0]/size[1]) + Math.sqrt(size[1]/size[0]) );
			Logger.printUserLog("Sample sizes meta-analysis 1: " + size[0]);
			Logger.printUserLog("Sample sizes meta-analysis 2: " + size[1]);
		}
		else
		{
			Logger.printUserLog("Summary statistics analysis for case-contrl studies.");
			double[] size = lamArgs.getCCsize();
			R1 = size[0]/size[1];
			R2 = size[2]/size[3];
			Logger.printUserLog("Sample size for meta-analysis 1: " + size[0] + " cases, " + size[1] + " controls. (R1=" + R1 + ")");
			Logger.printUserLog("Sample size for meta-analysis 2: " + size[2] + " cases, " + size[3] + " controls. (R2=" + R2 + ")");
			double s1 = size[0] + size[1];
			double s2 = size[2] + size[3];
			Kappa = 2 / (Math.sqrt(s1 /s2) + Math.sqrt(s2 / s1));
		}
		Logger.printUserLog("Kappa: " + Kappa);
	}

	private HashMap<String, MetaStat> readMeta(int metaIdx)
	{
		BufferedReader reader = BufferedReader.openTextFile(MetaFile[metaIdx], "Summary Statistic file");
		String[] tokens = reader.readTokens();
		int tokenLen = tokens.length;
		
		for(int i = 0; i < tokens.length; i++)
		{
			if (tokens[i].equalsIgnoreCase("snp"))
			{
				SMidx[metaIdx][0] = i;
			}
			if (lamArgs.isQT()) 
			{
				if (tokens[i].equalsIgnoreCase("beta"))
				{
					SMidx[metaIdx][1] = i;
				}
			}
			else
			{
				if (tokens[i].equalsIgnoreCase("or"))
				{
					SMidx[metaIdx][1] = i;
				}
			}
			if (tokens[i].equalsIgnoreCase("se"))
			{
				SMidx[metaIdx][2] = i;
			}
			if (tokens[i].equalsIgnoreCase("a1"))
			{
				SMidx[metaIdx][3] = i;
			}
			if (tokens[i].equalsIgnoreCase("a2"))
            {
				SMidx[metaIdx][4] = i;
			}
		}

		boolean qFlag = false;
		
		if (SMidx[metaIdx][0] == -1)
		{
			Logger.printUserLog("Cannot find the snp column in " + MetaFile[metaIdx]);
			qFlag = true;
		}
		if (SMidx[metaIdx][1] == -1)
		{
			Logger.printUserLog("Cannot find the effect column in " + MetaFile[metaIdx]);
			qFlag = true;
		}
		if (SMidx[metaIdx][2] == -1)
		{
			Logger.printUserLog("Cannot find the se column in " + MetaFile[metaIdx]);
			qFlag = true;
		}
		if (SMidx[metaIdx][3] == -1)
		{
			Logger.printUserLog("Cannot find the allele 1 in " + MetaFile[metaIdx]);
		}
//		if (SMidx[metaIdx][4] == -1)
//		{
//			Logger.printUserLog("Cannot find the allele 2 in " + MetaFile[metaIdx]);
//		}
		if (qFlag)
		{
			System.exit(0);
		}

		HashMap<String, MetaStat> sumstat = NewIt.newHashMap();
		int cnt = 0;
		int cntBadEffect = 0;
		int cntBadSE = 0;
		int cntBadA1 = 0;
		int cntBadA2 = 0;
		while( (tokens = reader.readTokens(tokenLen)) != null)
		{
			if (ConstValues.isNA(tokens[SMidx[metaIdx][1]]))
			{
				cntBadEffect++;
				continue;
			}
			if (ConstValues.isNA(tokens[SMidx[metaIdx][2]]))
			{
				cntBadSE++;
				continue;
			}
			if (Float.parseFloat(tokens[SMidx[metaIdx][2]]) == 0)
			{
				cntBadSE++;
				continue;
			}
			if (tokens[SMidx[metaIdx][3]].length() != 1)
			{
				cntBadA1++;
				continue;
			}
			if (SMidx[metaIdx][4] != -1)
			{
				if (tokens[SMidx[metaIdx][4]].length() != 1)
				{
					cntBadA2++;
					continue;
				}				
			}

			MetaStat ms = null;
			if (SMidx[metaIdx][4] == -1)
			{
				ms = new MetaStat(tokens[SMidx[metaIdx][0]], Float.parseFloat(tokens[SMidx[metaIdx][1]]), Float.parseFloat(tokens[SMidx[metaIdx][2]]), tokens[SMidx[metaIdx][3]].charAt(0));
			}
			else
			{
				ms = new MetaStat(tokens[SMidx[metaIdx][0]], Float.parseFloat(tokens[SMidx[metaIdx][1]]), Float.parseFloat(tokens[SMidx[metaIdx][2]]), tokens[SMidx[metaIdx][3]].charAt(0), tokens[SMidx[metaIdx][4]].charAt(0));				
			}
			sumstat.put(ms.getSNP(), ms);
			cnt++;
		}
		if(cntBadEffect > 0)
		{
			if (cntBadEffect == 1)
			{
				Logger.printUserLog("Removed " + cntBadEffect + " locus due to not numeric effect.");				
			}
			else
			{
				Logger.printUserLog("Removed " + cntBadEffect + " locus due to not numeric effect.");				
			}
		}
		if(cntBadSE > 0)
		{
			if (cntBadSE == 1)
			{
				Logger.printUserLog("Removed " + cntBadSE + " locus due to not numeric se.");				
			}
			else
			{
				Logger.printUserLog("Removed " + cntBadSE + " loci due to not numeric se.");
			}
		}
		if(cntBadA1 > 0)
		{
			if (cntBadA1 == 1)
			{
				Logger.printUserLog("Removed " + cntBadA1 + " locus due to bad a1 allele.");				
			}
			else
			{
				Logger.printUserLog("Removed " + cntBadA1 + " loci due to bad a1 allele.");
			}
		}
		if(cntBadA2 > 0)
		{
			if (cntBadA2 == 1)
			{
				Logger.printUserLog("Removed " + cntBadA2 + " locus due to bad a2 allele.");
			}
			else
			{
				Logger.printUserLog("Removed " + cntBadA2 + " loci due to bad a2 allele.");				
			}
		}
		if(cnt == 0)
		{
			Logger.printUserLog("Did not find any summary statistics from " + MetaFile[metaIdx]);
		}
		else
		{
			Logger.printUserLog("Read " + cnt + " summary statistics from " + MetaFile[metaIdx]);			
		}

		return sumstat;
	}

	private void calculateLambdaD()
	{
		ArrayList<Double> lD = NewIt.newArrayList();
		int cntAmbiguous = 0;
		for(Map.Entry<String, MetaStat> entry : SumStat1.entrySet())
		{
			String key = entry.getKey();
			if (!SumStat2.containsKey(key))
			{
				continue;
			}
			MetaStat ms1 = entry.getValue();
			MetaStat ms2 = SumStat2.get(key);
			double d = 0;
			
			if (SMidx[0][4] != -1)
			{
				SNPMatch.isAmbiguous(ms1.getA1(), ms1.getA2());
				cntAmbiguous++;
				continue;
			}
			if (SMidx[1][4] != -1)
			{
				SNPMatch.isAmbiguous(ms1.getA1(), ms2.getA2());
				cntAmbiguous++;
				continue;
			}
			
			if (ms1.getA1() == ms2.getA1() || ms1.getA1() == SNPMatch.Flip(ms2.getA1()))
			{
				d = (ms1.getEffect() - ms2.getEffect()) * (ms1.getEffect() - ms2.getEffect()) / (ms1.getSE() * ms1.getSE() + ms2.getSE() * ms2.getSE());
			}
			else
			{
				if (lamArgs.isQT())
				{
					d = ((-1) * ms1.getEffect() - ms2.getEffect()) * ((-1) * ms1.getEffect() - ms2.getEffect()) / (ms1.getSE() * ms1.getSE() + ms2.getSE() * ms2.getSE());
				}
				else
				{
					d = (1/ ms1.getEffect() - ms2.getEffect()) * (1 / ms1.getEffect() - ms2.getEffect()) / (ms1.getSE() * ms1.getSE() + ms2.getSE() * ms2.getSE());					
				}
			}
			lD.add(d);
        }
		if (cntAmbiguous > 0)
		{
			if (cntAmbiguous == 1)
			{
				Logger.printUserLog("Removed " + cntAmbiguous + " locus (AT/GC).");				
			}
			else
			{
				Logger.printUserLog("Removed " + cntAmbiguous + " loci (AT/GC).");				
			}
		}
		Logger.printUserLog("Lambda is calculated based on " + (lD.size() - cntAmbiguous) + " summary statistics between two files.");
		double[] ld = new double[lD.size()];
		for(int i = 0; i < ld.length; i++)
		{
			ld[i] = lD.get(i).doubleValue();
		}

		DescriptiveStatistics ds = new DescriptiveStatistics(ld);
		double[] Sortld = ds.getSortedValues();
		if (ds.getN() % 2 == 0) 
		{
			LambdaMedian = (Sortld[(int) (ds.getN()/2)] + Sortld[(int) (ds.getN()/2) - 1])/2 / 0.4549364;
		}
		else
		{
			LambdaMedian = Sortld[(int) ((ds.getN()-1)/2)] / 0.4549364;
		}
		LambdaMean = ds.getMean();
		rhoMedian = (1-LambdaMedian) / Kappa;
		rhoMean = (1 - LambdaMean) / Kappa;

		if (lamArgs.isQT())
		{
			double[] qtSize = lamArgs.getQTsize();
			OSMedian = (1 - LambdaMedian) / Kappa * Math.sqrt(qtSize[0] * qtSize[1]);
			OSMean = (1 - LambdaMean) / Kappa * Math.sqrt(qtSize[0] * qtSize[1]);
		}
		else
		{
			double[] ccSize = lamArgs.getCCsize();
			OSMedian = (1 - LambdaMedian) / Kappa * Math.sqrt( (ccSize[0] + ccSize[1] ) * (ccSize[2] + ccSize[3]));
			OSMean = (1 - LambdaMean) / Kappa * Math.sqrt((ccSize[0] + ccSize[1] ) * (ccSize[2] + ccSize[3]));
			OSCtrlMedian = OSMedian / Math.sqrt(R1 * R2);
			OSCtrlMean = OSMean / Math.sqrt(R1 * R2);
		}
	}

	private double R1 = 1;
	private double R2 = 1;
	private double Kappa = 1;
	private LambdaDCommandArguments lamArgs;

	private int[][] SMidx = { {-1,-1,-1,-1,-1}, {-1,-1,-1,-1,-1}}; //snp, beta, se, a1, a2
	private String[] MetaFile = {null, null};
	private HashMap<String, MetaStat> SumStat1;
	private HashMap<String, MetaStat> SumStat2;

	private double LambdaMedian = 0;
	private double LambdaMean = 0;
	private double OSMedian = 0;
	private double OSMean = 0;
	private double OSCtrlMedian = 0;
	private double OSCtrlMean = 0;
	private double rhoMedian = 0;
	private double rhoMean = 0;
}
