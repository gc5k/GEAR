package gear.subcommands.exsnp;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;

import gear.family.pedigree.file.SNP;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.BufferedReader;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;

public class ExSNPCommandImpl  extends CommandImpl
{

	@Override
	public void execute(CommandArguments cmdArgs)
	{
		ExSNPCommandArguments esArgs = (ExSNPCommandArguments) cmdArgs;
		ArrayList<String> bFile = esArgs.getBFiles();
		if(bFile.size() == 0)
		{
			Logger.printUserLog("no bfiles found. GEAR quitted.");
			System.exit(0);
		}
		for(int i = 0; i < bFile.size(); i++)
		{
			int SNPcnt = 0;
			String bim = bFile.get(i) + ".bim";
			BufferedReader reader = BufferedReader.openTextFile(bim, "Extract SNP Batch");
			String[] tokens = null;
			while((tokens = reader.readTokens())!=null)
			{
				if(tokens.length == 6)
				{
					if (i == 0)
					{
						SNP snp = new SNP(tokens[0], tokens[1], Float.parseFloat(tokens[2]), Integer.parseInt(tokens[3]), tokens[4].charAt(0), tokens[5].charAt(0));
						snpInfo.add(snp);
					}
					if (snpCnt.containsKey(tokens[1]))
					{
						Integer cnt = snpCnt.get(tokens[1]);
						cnt++;
						snpCnt.put(tokens[1], cnt);
					}
					else
					{
						snpCnt.put(tokens[1], 1);
					}
				}
				SNPcnt++;
			}
			Logger.printUserLog("Read " + SNPcnt + " snps in '" + bFile.get(i) + "'.");
		}

		ArrayList<SNP> finalSNP = NewIt.newArrayList();

		for(int i = 0; i < snpInfo.size(); i++)
		{
			SNP snp = snpInfo.get(i);
			String key = snp.getName();
			Integer cnt = snpCnt.get(key);
			if(cnt.intValue() == bFile.size())
			{
				finalSNP.add(snp);
			}
		}

		Collections.sort(finalSNP);
		
		StringBuffer sb = new StringBuffer();
		sb.append(esArgs.getOutRoot());
		sb.append(".comsnp");
		PrintStream ps = FileUtil.CreatePrintStream(sb.toString());

		for(int i = 0; i < finalSNP.size(); i++)
		{
			ps.println(finalSNP.get(i).getName());
		}
		ps.close();
	}

	private HashMap<String, Integer> snpCnt = NewIt.newHashMap();
	private ArrayList<SNP> snpInfo = NewIt.newArrayList();
}
