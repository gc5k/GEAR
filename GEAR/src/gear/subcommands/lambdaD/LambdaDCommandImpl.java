package gear.subcommands.lambdaD;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import gear.ConstValues;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.BufferedReader;
import gear.util.Logger;
import gear.util.NewIt;
import gear.util.SNPMatch;

public class LambdaDCommandImpl extends CommandImpl
{
	private void initial()
	{
		MetaFile = lamArgs.getMetaFile();
		mat = new double[MetaFile.length][MetaFile.length];
		lamMat = new double[MetaFile.length][MetaFile.length];
		olMat = new double[MetaFile.length][MetaFile.length];
		olCtrlMat = new double[MetaFile.length][MetaFile.length];
		kMat = new double[MetaFile.length][MetaFile.length];

		if (MetaFile.length < 2)
		{
			Logger.printUserError("At least two meta files should be specified.");
			Logger.printUserError("GEAR quitted.");
			System.exit(0);
		}

		logit = new boolean[MetaFile.length];
		Arrays.fill(logit, false);

		KeyIdx = new int[MetaFile.length][5];
		for (int i = 0; i < KeyIdx.length; i++)
		{
			Arrays.fill(KeyIdx[i], -1);
		}
		
		for (int i = 0; i < MetaFile.length; i++)
		{
			HashMap<String, MetaStat> m = readMeta(i);
			meta.add(m);
		}
	}

	@Override
	public void execute(CommandArguments cmdArgs)
	{
		lamArgs = (LambdaDCommandArguments) cmdArgs;
		initial();

		for (int i=0; i < MetaFile.length -1; i++)
		{
			mat[i][i] = 1;
//			SumStat1 = readMeta(i);
			for (int j = (i+1); j < MetaFile.length; j++)
			{
//				SumStat2 = readMeta(j);

				if (lamArgs.isQT())
				{
					Logger.printUserLog("Summary statistics analysis for quantitative traits.");
					double[] size = lamArgs.getQTsize();
					Kappa = 2 / ( Math.sqrt(size[i]/size[j]) + Math.sqrt(size[j]/size[i]) );
					Logger.printUserLog("Sample sizes for '" + MetaFile[i] + "': " + size[0]);
					Logger.printUserLog("Sample sizes for '" + MetaFile[j] + "': " +size[1]);
				}
				else
				{
					Logger.printUserLog("Summary statistics analysis for case-contrl studies.");
					double[] size = lamArgs.getCCsize();
					R1 = size[i*2]/size[i*2+1];
					R2 = size[j*2]/size[j*2+1];
					Logger.printUserLog("Sample size for '" + MetaFile[i] + "': " + size[i*2] + " cases, " + size[i*2+1] + " controls; R1 = " + R1 + ".");
					Logger.printUserLog("Sample size for '" + MetaFile[j] + "': " + size[j*2] + " cases, " + size[j*2+1] + " controls; R2 = " + R2 + ".");
					double s1 = size[i*2] + size[i*2+1];
					double s2 = size[j*2] + size[j*2+1];
					Kappa = 2 / (Math.sqrt(s1 / s2) + Math.sqrt(s2 / s1));
				}
				Logger.printUserLog("Kappa: " + Kappa);

				calculateLambdaD(i, j);
				Logger.printUserLog("LambdaD median: " + LambdaMedian);
				Logger.printUserLog("Estimated rho (lambdaD median): " + rhoMedian);
				Logger.printUserLog("Estimated overlapping samples (lambdaD median): " + OSMedian);
				if (!lamArgs.isQT())
				{
					Logger.printUserLog("Estiamted overlapping controls: " + OSCtrlMedian);
				}
				Logger.printUserLog("LambdaD mean: " + LambdaMean);
				Logger.printUserLog("Estimated rho (lambdaD mean): " + rhoMean);
				Logger.printUserLog("Estimated overlapping samples (lambdaD mean): " + OSMean);
				if (!lamArgs.isQT())
				{
					Logger.printUserLog("Estiamted overlapping controls: " + OSCtrlMean);
				}
				Logger.printUserLog("\n");

				mat[i][j] = mat[j][i] = rhoMedian;
				lamMat[i][j] = lamMat[j][i] = LambdaMedian;
				olMat[i][j] = olMat[j][i] = OSMedian;
				kMat[i][j] = kMat[j][i] = Kappa;
				if (!lamArgs.isQT())
				{
					olCtrlMat[i][j] = olCtrlMat[j][i] = OSCtrlMedian;
				}
//				printOut(i,j);
			}
		}

		WriteMat();
	}

	private HashMap<String, MetaStat> readMeta(int metaIdx)
	{
		BufferedReader reader = null;
		if (lamArgs.isGZ())
		{
			reader = BufferedReader.openGZipFile(MetaFile[metaIdx], "Summary Statistic file");
			
		}
		else
		{
			reader = BufferedReader.openTextFile(MetaFile[metaIdx], "Summary Statistic file");
		}

		String[] tokens = reader.readTokens();
		int tokenLen = tokens.length;

		for(int i = 0; i < tokens.length; i++)
		{
			if (tokens[i].equalsIgnoreCase(lamArgs.getKey(LambdaDCommandArguments.SNP)))
			{
				KeyIdx[metaIdx][0] = i;
			}
			if (lamArgs.isQT()) 
			{
				if (tokens[i].equalsIgnoreCase(lamArgs.getKey(LambdaDCommandArguments.BETA)))
				{
					KeyIdx[metaIdx][1] = i;
				}
			}
			else
			{
				if(tokens[i].equalsIgnoreCase(lamArgs.getKey(LambdaDCommandArguments.BETA)))
				{
					KeyIdx[metaIdx][1] = i;
					logit[metaIdx] = false;
				}
				else if (tokens[i].equalsIgnoreCase(lamArgs.getKey(LambdaDCommandArguments.OR)))
				{
					KeyIdx[metaIdx][1] = i;
					logit[metaIdx] = true;
				} 
			}
			if (tokens[i].equalsIgnoreCase(lamArgs.getKey(LambdaDCommandArguments.SE)))
			{
				KeyIdx[metaIdx][2] = i;
			}
			if (tokens[i].equalsIgnoreCase(lamArgs.getKey(LambdaDCommandArguments.A1)))
			{
				KeyIdx[metaIdx][3] = i;
			}
			if (tokens[i].equalsIgnoreCase(lamArgs.getKey(LambdaDCommandArguments.A2)))
            {
				KeyIdx[metaIdx][4] = i;
			}
		}

		boolean qFlag = false;

		if (KeyIdx[metaIdx][0] == -1)
		{
			Logger.printUserLog("Cannot find the snp column in " + MetaFile[metaIdx]);
			qFlag = true;
		}
		if (KeyIdx[metaIdx][1] == -1)
		{
			Logger.printUserLog("Cannot find the effect column in " + MetaFile[metaIdx]);
			qFlag = true;
		}
		if (KeyIdx[metaIdx][2] == -1)
		{
			Logger.printUserLog("Cannot find the se column in " + MetaFile[metaIdx]);
			qFlag = true;
		}
		if (KeyIdx[metaIdx][3] == -1)
		{
			Logger.printUserLog("Cannot find the allele 1 in " + MetaFile[metaIdx]);
		}
//		if (SMidx[metaIdx][4] == -1)
//		{
//			Logger.printUserLog("Cannot find the allele 2 in " + MetaFile[metaIdx]);
//		}
		if (qFlag)
		{
			Logger.printUserLog("GEAR quitted.");
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
			if (ConstValues.isNA(tokens[KeyIdx[metaIdx][1]]))
			{
				cntBadEffect++;
				continue;
			}
			if (ConstValues.isNA(tokens[KeyIdx[metaIdx][2]]))
			{
				cntBadSE++;
				continue;
			}
			if (Float.parseFloat(tokens[KeyIdx[metaIdx][2]]) == 0)
			{
				cntBadSE++;
				continue;
			}
			if (tokens[KeyIdx[metaIdx][3]].length() != 1)
			{
				cntBadA1++;
				continue;
			}
			if (KeyIdx[metaIdx][4] != -1)
			{
				if (tokens[KeyIdx[metaIdx][4]].length() != 1)
				{
					cntBadA2++;
					continue;
				}
			}

			MetaStat ms = null;
			if (KeyIdx[metaIdx][4] == -1)
			{
				ms = new MetaStat(tokens[KeyIdx[metaIdx][0]], Float.parseFloat(tokens[KeyIdx[metaIdx][1]]), Float.parseFloat(tokens[KeyIdx[metaIdx][2]]), tokens[KeyIdx[metaIdx][3]].charAt(0), logit[metaIdx]);
			}
			else
			{
				ms = new MetaStat(tokens[KeyIdx[metaIdx][0]], Float.parseFloat(tokens[KeyIdx[metaIdx][1]]), Float.parseFloat(tokens[KeyIdx[metaIdx][2]]), tokens[KeyIdx[metaIdx][3]].charAt(0), tokens[KeyIdx[metaIdx][4]].charAt(0), logit[metaIdx]);			
			}
			sumstat.put(ms.getSNP(), ms);
			cnt++;
		}
		if(cnt == 0)
		{
			Logger.printUserLog("Did not find any summary statistics from '" + MetaFile[metaIdx] + "'");
		}
		else
		{
			Logger.printUserLog("Read " + cnt + " effective summary statistics from '" + MetaFile[metaIdx] + "'");			
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

		return sumstat;
	}

	private void calculateLambdaD(int idx1, int idx2)
	{
		ArrayList<Double> lD = NewIt.newArrayList();
		int cntAmbiguous = 0;
		HashMap<String, MetaStat> SumStat1 = meta.get(idx1);
		HashMap<String, MetaStat> SumStat2 = meta.get(idx2);

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

			if (KeyIdx[idx1][4] != -1)
			{
				if (SNPMatch.isAmbiguous(ms1.getA1(), ms1.getA2()))
				{
					cntAmbiguous++;
					continue;
				}
			}
			if (KeyIdx[idx2][4] != -1)
			{
				if (SNPMatch.isAmbiguous(ms2.getA1(), ms2.getA2()))
				{
					cntAmbiguous++;
					continue;					
				}
			}

			if (ms1.getA1() == ms2.getA1() || ms1.getA1() == SNPMatch.Flip(ms2.getA1()))
			{
				d = (ms1.getEffect() - ms2.getEffect()) * (ms1.getEffect() - ms2.getEffect()) / (ms1.getSE() * ms1.getSE() + ms2.getSE() * ms2.getSE());
			}
			else
			{
				d = ((-1) * ms1.getEffect() - ms2.getEffect()) * ((-1) * ms1.getEffect() - ms2.getEffect()) / (ms1.getSE() * ms1.getSE() + ms2.getSE() * ms2.getSE());
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
		double mean = 0;
		for(int i = 0; i < ld.length; i++)
		{
			ld[i] = lD.get(i).doubleValue();
			mean += ld[i];
		}

		Arrays.sort(ld);
		if (ld.length % 2 == 0)
		{
			LambdaMedian = (ld[(int) (ld.length/2)] + ld[(int) (ld.length/2) - 1])/2 / 0.4549364;
		}
		else
		{
			LambdaMedian = ld[(int) ((ld.length-1)/2)] / 0.4549364;
		}

		LambdaMean = mean/ld.length;
		rhoMedian = (1-LambdaMedian) / Kappa;
		rhoMean = (1 - LambdaMean) / Kappa;

		if (lamArgs.isQT())
		{
			double[] qtSize = lamArgs.getQTsize();
			OSMedian = (1 - LambdaMedian) / Kappa * Math.sqrt(qtSize[idx1] * qtSize[idx2]);
			OSMean = (1 - LambdaMean) / Kappa * Math.sqrt(qtSize[idx1] * qtSize[idx2]);
		}
		else
		{
			double[] ccSize = lamArgs.getCCsize();
			OSMedian = (1 - LambdaMedian) / Kappa * Math.sqrt( (ccSize[idx1*2] + ccSize[idx1*2+1] ) * (ccSize[idx2*2] + ccSize[idx2*2+1]) );
			OSMean = (1 - LambdaMean) / Kappa * Math.sqrt( (ccSize[idx1*2] + ccSize[idx1*2+1] ) * (ccSize[idx2*2] + ccSize[idx2*2+1]) );
			OSCtrlMedian = OSMedian / Math.sqrt(R1 * R2);
			OSCtrlMean = OSMean / Math.sqrt(R1 * R2);
		}
	}

	private void printOut(int idx1, int idx2)
	{
		PrintWriter writer = null;
		try 
		{
			writer = new PrintWriter(new BufferedWriter(new FileWriter(lamArgs.getOutRoot() + "." + (idx1+1) + "-" + (idx2+1) + ".lam")));
		}
		catch (IOException e)
		{
			Logger.handleException(e, "An I/O exception occurred when writing '" + lamArgs.getOutRoot() + ".lam" + "'.");
		}
		writer.println("Meta File " + (idx1+1) + ": " +MetaFile[idx1]);
		writer.println("Meta File " + (idx2+1) + ": " +MetaFile[idx2]);

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
			writer.println("Overlapping controls: " + OSCtrlMedian);
		}
		writer.close();
	}

	private void WriteMat()
	{
		PrintWriter writer = null;
		try 
		{
			writer = new PrintWriter(new BufferedWriter(new FileWriter(lamArgs.getOutRoot() + ".lmat")));
		}
		catch (IOException e)
		{
			Logger.handleException(e, "An I/O exception occurred when writing '" + lamArgs.getOutRoot() + ".lmat" + "'.");
		}
		
		writer.println("correlation:");
		for (int i = 0; i < mat.length; i++)
		{
			for (int j = 0; j < mat[i].length; j++)
			{
				writer.print(String.format("%.4f", mat[i][j]) + " ");
			}
			writer.println();
		}

		writer.println("LambdaD:");
		for (int i = 0; i < lamMat.length; i++)
		{
			for (int j = 0; j < lamMat[i].length; j++)
			{
				writer.print(String.format("%.4f", lamMat[i][j]) + " ");
			}
			writer.println();
		}

		writer.println("Kappa:");
		for (int i = 0; i < kMat.length; i++)
		{
			for (int j = 0; j < kMat[i].length; j++)
			{
				writer.print(String.format("%.4f", kMat[i][j]) + " ");
			}
			writer.println();
		}

		writer.println("Overlapping Samples:");
		for (int i = 0; i < olMat.length; i++)
		{
			for (int j = 0; j < olMat[i].length; j++)
			{
				writer.print(String.format("%.4f", olMat[i][j]) + " ");
			}
			writer.println();
		}
		if(!lamArgs.isQT())
		{
			writer.println("Overlapping Controls:");
			for (int i = 0; i < olCtrlMat.length; i++)
			{
				for (int j = 0; j < olCtrlMat[i].length; j++)
				{
					writer.print(String.format("%.4f", olCtrlMat[i][j]) + " ");
				}
				writer.println();
			}
		}
		
		writer.close();
	}

	private double R1 = 1;
	private double R2 = 1;
	private double Kappa = 1;
	private LambdaDCommandArguments lamArgs;

	private int[][] KeyIdx; //snp, beta, se, a1, a2
	private String[] MetaFile;
	private ArrayList<HashMap<String, MetaStat>> meta = NewIt.newArrayList();
//	private HashMap<String, MetaStat> SumStat1;
//	private HashMap<String, MetaStat> SumStat2;

	private double LambdaMedian = 0;
	private double LambdaMean = 0;
	private double OSMedian = 0;
	private double OSMean = 0;
	private double OSCtrlMedian = 0;
	private double OSCtrlMean = 0;
	private double rhoMedian = 0;
	private double rhoMean = 0;

	private boolean[] logit;
	
	private double[][] mat;
	private double[][] lamMat;
	private double[][] olMat;
	private double[][] olCtrlMat;
	private double[][] kMat;
}
