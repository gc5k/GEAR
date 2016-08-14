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
import gear.subcommands.weightedmeta.util.GMRes;
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
		readNSS();
		
		oath();
	}
	
	private void oath()
	{
		Logger.printUserLog("Starting OATH-analysis...");
		int totalCnt = 0;
		int cnt = 0;
		int singularCnt = 0;
		int atgcCnt = 0;
		
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
		Logger.printUserLog("In total "+ cnt + " loci have been used for meta-analysis.");
		if (singularCnt > 0)
		{
			Logger.printUserLog(singularCnt + " loci were excluded from analyais because of singular matrix.");
		}
		PrintSynResults();
	}

	private void readCM()
	{
		Logger.printUserLog("Reading correlation matrix from '" + synArgs.getCMFile() + "'.");
		BufferedReader reader = null;
		reader = BufferedReader.openTextFile(synArgs.getCMFile(), "Correlation matrix file.");

		String[] tokens = reader.readTokens();
		int tokenLen = tokens.length;

		if (tokenLen != synArgs.getMetaFile().length)
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
	}

	private void readNSS() 
	{
		boolean[] FileKeep = new boolean[synArgs.getMetaFile().length];
		Arrays.fill(FileKeep, true);
		fReader = new SynthFReader(synArgs.getMetaFile(), FileKeep,
				synArgs.getKeys(), true, synArgs.isGZ(),
				synArgs.isChr(), synArgs.getChr());

		fReader.Start();

		int NumMetaFile = synArgs.getMetaFile().length;

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
			Logger.handleException(e, "An I/O exception occurred when writing '" + synArgs.getOutRoot() + ".gmeta" + "'.\n");
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

	private SynthCommandArguments synArgs;
	private SynthFReader fReader;
	private double[][] corMat;
	private ArrayList<SynthRes> grArray = NewIt.newArrayList();

}
