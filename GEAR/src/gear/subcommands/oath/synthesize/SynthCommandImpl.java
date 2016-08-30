package gear.subcommands.oath.synthesize;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.Set;

import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.subcommands.oath.synthesize.freader.SynthFReader;
import gear.subcommands.oath.synthesize.freader.SynthFStat;
import gear.subcommands.oath.synthesize.freader.SynthMatrix;
import gear.subcommands.oath.synthesize.freader.SynthRes;
import gear.util.BufferedReader;
import gear.util.Logger;
import gear.util.NewIt;

public class SynthCommandImpl extends CommandImpl 
{

	@Override
	public void execute(CommandArguments cmdArgs) 
	{
		synArgs = (SynthCommandArguments) cmdArgs;
		readCM();
		readKeep();
		if (synArgs.getKeepBatchIdx() == null && corMat.length != synArgs.getNSSFile().length)
		{
			Logger.printUserLog("The size of correlation matrix ("+ corMat.length +"X"+corMat.length+")does not fit to the nss files.");
			System.exit(1);
		}
		readNSS();
		oath();
	}
	
	private void readKeep()
	{
		keepBatch = new boolean[corMat.length];
		int[] kp = synArgs.getKeepBatchIdx();
		if (kp != null)
		{
			Arrays.sort(kp);
			Arrays.fill(keepBatch, false);
			for(int i = 0; i < kp.length; i++)
			{
				if (kp[i] > keepBatch.length) 
				{
					Logger.printUserLog("--keep-nss index is greater than " + keepBatch.length + ".");
					Logger.printUserLog("GEAR quitted.");
					System.exit(1);
				}
				keepBatch[kp[i]] = true;
			}

			double[][] cm = new double[kp.length][kp.length];
			for(int i = 0; i < kp.length; i++)
			{
				for(int j = 0; j < kp.length; j++)
				{
					cm[i][j] = corMat[kp[i]][kp[j]];
				}
			}
			corMat = cm;
			Logger.printUserLog(kp.length + " nss files are kept for analysis.");
		}
		else
		{
			Arrays.fill(keepBatch, true);
		}
	}

	private void oath()
	{
		Logger.printUserLog("");
		Logger.printUserLog("Starting OATH-analysis...");
		int totalCnt = 0;
		int singularCnt = 0;
//		int atgcCnt = 0;
		
		Set<String> snps = fReader.getMetaSNPTable().keySet();
		for (Iterator<String> e=snps.iterator(); e.hasNext();)
		{
			String snp = e.next();
//			System.out.println(snp);
			ArrayList<Integer> Int = fReader.getMetaSNPTable().get(snp);

			SynthFStat ms = null;
//			int i = 0;
			for(int i = 0; i < (Int.size() - 1); i++)
			{
				if(Int.get(i).intValue() != 0) break; 
			}
			ms = fReader.getMetaStat().get(0).get(snp);

			if (Int.get(Int.size()-1).intValue() != (Int.size() -1))
			{
				continue;
			}
			SynthMatrix SynMat = new SynthMatrix(synArgs.getN(), snp, Int, corMat, fReader);

			SynthRes sRes = new SynthRes(Int.get(Int.size() - 1).intValue());
			sRes.SetSNP(snp);
			sRes.SetChr(ms.getChr());
			sRes.SetBP(ms.getBP());
			sRes.SetA1(ms.getA1());
			sRes.SetA2(ms.getA2());

			sRes.SetB(SynMat.getOATHB());
			sRes.SetSE(SynMat.getOATHse());
			sRes.SetZ(SynMat.getOATHB()/SynMat.getOATHse());
			totalCnt++;
			grArray.add(sRes);
		}
		Collections.sort(grArray);

		Logger.printUserLog("In total "+ totalCnt + " loci have been read.");
		Logger.printUserLog("In total "+ grArray.size() + " loci have been used for OATH analysis.");
		if (singularCnt > 0)
		{
			Logger.printUserLog(singularCnt + " loci were excluded from analyais because of singular matrix.");
		}
		PrintSynResults();
	}

	private void readCM()
	{
		BufferedReader reader = null;
		reader = BufferedReader.openTextFile(synArgs.getCMFile(), "Correlation matrix file.");

		String[] tokens = reader.readTokens();
		int tokenLen = tokens.length;

		if (tokenLen != synArgs.getNSSFile().length)
		{
			Logger.printUserError("The dimension of the matrix does not match the number of nss files.");
			Logger.printUserError("Gear quitted.");
			System.exit(1);
		}
		corMat = new double[tokenLen][tokenLen];

		int cnt=0;
		do
		{
			for(int i = 0; i < tokens.length; i++)
			{
				corMat[cnt][i] = Double.parseDouble(tokens[i]);
			}
			cnt++;
		} while ((tokens = reader.readTokens()) != null);
		Logger.printUserLog("Reading " +tokenLen +"X" +tokenLen +" correlation matrix from '" + synArgs.getCMFile() + "'.");
	}

	private void readNSS() 
	{
		fReader = new SynthFReader(synArgs.getNSSFile(), keepBatch,
				synArgs.getKeys(), true, synArgs.isGZ(),
				synArgs.isChr(), synArgs.getChr());

		fReader.Start();

		int NumMetaFile = synArgs.getNSSFile().length;

		if (NumMetaFile < 2)
		{
			Logger.printUserError("At least two summary statistic files should be specified.\n");
			Logger.printUserError("GEAR quitted.\n");
			System.exit(0);
		}
	}

	private void PrintSynResults()
	{
        PrintWriter writer = null;
        try
        {
        	writer = new PrintWriter(new BufferedWriter(new FileWriter(synArgs.getOutRoot()+".oath")));
        	Logger.printUserLog("Writting detailed test statistics into '"+synArgs.getOutRoot() + ".oath.'\n");
        }
		catch (IOException e)
		{
			Logger.handleException(e, "An I/O exception occurred when writing '" + synArgs.getOutRoot() + ".oath'.\n");
		}

		for(int i = 0; i < grArray.size(); i++)
		{
			SynthRes gr = grArray.get(i);
			if(i == 0)
			{
				writer.write(gr.printTitle() + "\n");
			}
			writer.write(gr.toString()+ "\n");
		}
		writer.close();
	}

	private boolean[] keepBatch = null;
	private SynthCommandArguments synArgs;
	private SynthFReader fReader;
	private double[][] corMat = null;
	private ArrayList<SynthRes> grArray = NewIt.newArrayList();

}
