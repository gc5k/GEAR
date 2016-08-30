package gear.subcommands.oath.synthesize.freader;

import gear.ConstValues;
import gear.util.BufferedReader;
import gear.util.Logger;
import gear.util.NewIt;
import gear.subcommands.oath.OATHConst;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

public class SynthFReader
{

	public SynthFReader(String[] MetaFile, boolean[] FileKeep, String[] field, boolean isQT, boolean isGZ, boolean isChr, int Chr)
	{
		this.field = field;
		this.isGZ = isGZ;

		workingMetaFile = NewIt.newArrayList();
		this.isChr = isChr;
		this.chrKeep = Chr;

		for (int i = 0; i < FileKeep.length; i++)
		{
			if (FileKeep[i])
			{
				workingMetaFile.add(MetaFile[i]);
			}
		}

		if (workingMetaFile.size() == 0)
		{
			Logger.printUserLog("No cohort left. GEAR quitted.");
			System.exit(0);
		}
		else
		{
			Logger.printUserLog(workingMetaFile.size() + " cohorts are remained for analysis.");
		}

		KeyIdx = new int[workingMetaFile.size()][OATHConst.num_key];
		for (int i = 0; i < KeyIdx.length; i++)
		{
			Arrays.fill (KeyIdx[i], -1);
		}
	}

	public void Start()
	{
		for (int i = 0; i < workingMetaFile.size(); i++)
		{
			HashMap<String, SynthFStat> m = readMeta(i);
			MStat.add(m);
		}	
	}

	public int getCohortNum()
	{
		return workingMetaFile.size();
	}
	
	public String[] getMetaFile()
	{
		return workingMetaFile.toArray(new String[0]);
	}

	public int[][] getKeyIndex()
	{
		return KeyIdx;
	}

	public int getNumMetaFile()
	{
		return workingMetaFile.size();
	}

	public ArrayList<HashMap<String, SynthFStat>> getMetaStat()
	{
		return MStat;
	}

	public ArrayList<ArrayList<String>> getMetaSNPArray() 
	{
		return MetaSNPArray;
	}

	public HashMap<String, ArrayList<Integer>> getMetaSNPTable()
	{
		return MetaSNPTable;
	}

	private HashMap<String, SynthFStat> readMeta(int metaIdx)
	{
		BufferedReader reader = null;
		if (isGZ)
		{
			reader = BufferedReader.openGZipFile(workingMetaFile.get(metaIdx), "Summary Statistic file");
		}
		else
		{
			reader = BufferedReader.openTextFile(workingMetaFile.get(metaIdx), "Summary Statistic file");
		}

		String[] tokens = reader.readTokens();
		int tokenLen = tokens.length;

		for(int i = 0; i < tokens.length; i++)
		{
			if (tokens[i].equalsIgnoreCase(field[OATHConst.snp]))
			{
				KeyIdx[metaIdx][OATHConst.snp] = i;
			}
			if (tokens[i].equalsIgnoreCase(field[OATHConst.chr]))
			{
				KeyIdx[metaIdx][OATHConst.chr] = i;
			}
			if (tokens[i].equalsIgnoreCase(field[OATHConst.bp]))
			{
				KeyIdx[metaIdx][OATHConst.bp] = i;
			}
			if (tokens[i].equalsIgnoreCase(field[OATHConst.refale]))
			{
				KeyIdx[metaIdx][OATHConst.refale] = i;
			}
			if (tokens[i].equalsIgnoreCase(field[OATHConst.altale]))
            {
				KeyIdx[metaIdx][OATHConst.altale] = i;
			}
			if (tokens[i].equalsIgnoreCase(field[OATHConst.raf]))
			{
				KeyIdx[metaIdx][OATHConst.raf] = i;
			}
			if (tokens[i].equalsIgnoreCase(field[OATHConst.vg]))
			{
				KeyIdx[metaIdx][OATHConst.vg] = i;
			}
			if (tokens[i].equalsIgnoreCase(field[OATHConst.beta]))
			{
				KeyIdx[metaIdx][OATHConst.beta] = i;
			}
			if (tokens[i].equalsIgnoreCase(field[OATHConst.se]))
			{
				KeyIdx[metaIdx][OATHConst.se] = i;
			}
			if (tokens[i].equalsIgnoreCase(field[OATHConst.chi]))
			{
				KeyIdx[metaIdx][OATHConst.chi] = i;
			}
			if (tokens[i].equalsIgnoreCase(field[OATHConst.p]))
			{
				KeyIdx[metaIdx][OATHConst.p] = i;
			}
		}

		boolean qFlag = false;

		if (KeyIdx[metaIdx][OATHConst.snp] == -1)
		{
			Logger.printUserLog("Cannot find the snp column in " + workingMetaFile.get(metaIdx));
			qFlag = true;
		}

		if (KeyIdx[metaIdx][OATHConst.chr] == -1)
		{
			Logger.printUserLog("Cannot find the chr column in " + workingMetaFile.get(metaIdx));
			qFlag = true;
		}

		if (KeyIdx[metaIdx][OATHConst.bp] == -1)
		{
			Logger.printUserLog("Cannot find the bp column in " + workingMetaFile.get(metaIdx));
		}

		if (KeyIdx[metaIdx][OATHConst.refale] == -1)
		{
			Logger.printUserLog("Cannot find the allele 1 (reference) column in " + workingMetaFile.get(metaIdx));
			qFlag = true;
		}

		if (KeyIdx[metaIdx][OATHConst.altale] == -1)
		{
			Logger.printUserLog("Cannot find the allele 2 (alternative) column in " + workingMetaFile.get(metaIdx));
			qFlag = true;
		}

		if (KeyIdx[metaIdx][OATHConst.altale] == -1)
		{
			Logger.printUserLog("Cannot find the allele 2 (alternative) column in " + workingMetaFile.get(metaIdx));
			qFlag = true;
		}

		if (KeyIdx[metaIdx][OATHConst.raf] == -1)
		{
			Logger.printUserLog("Cannot find the allele frequency column in " + workingMetaFile.get(metaIdx));
			qFlag = true;
		}

		if (KeyIdx[metaIdx][OATHConst.vg] == -1)
		{
			Logger.printUserLog("Cannot find the allele frequency variance column in " + workingMetaFile.get(metaIdx));
			qFlag = true;
		}

		if (KeyIdx[metaIdx][OATHConst.beta] == -1)
		{
			Logger.printUserLog("Cannot find the beta column in " + workingMetaFile.get(metaIdx));
			qFlag = true;
		}

		if (KeyIdx[metaIdx][OATHConst.se] == -1)
		{
			Logger.printUserLog("Cannot find the se column in " + workingMetaFile.get(metaIdx));
			qFlag = true;
		}

		if (qFlag)
		{
			Logger.printUserLog("GEAR quitted.");
			System.exit(0);
		}

		HashMap<String, SynthFStat> sumstat = NewIt.newHashMap();
		ArrayList<String> snpArray = NewIt.newArrayList();
		int total = 0;
		int cnt = 0;
		int cntDup = 0;
		int cntBadChr = 0;
//		int cntBadBp = 0;
		int cntBadA1 = 0;
		int cntBadA2 = 0;
		int cntMissSNP = 0;
		int cntBadFreq = 0;
		int cntBadVg = 0;
		int cntBadBeta = 0;
		int cntBadSe = 0;

		while( (tokens = reader.readTokens(tokenLen)) != null)
		{
			total++;
			if (OATHConst.isNASNP(tokens[KeyIdx[metaIdx][OATHConst.snp]]))
			{
				cntMissSNP++;
				continue;
			}

			if (OATHConst.isNA(tokens[KeyIdx[metaIdx][OATHConst.chr]]))
			{
				cntBadChr++;
				continue;
			}

			if (OATHConst.isNA(tokens[KeyIdx[metaIdx][OATHConst.beta]]))
			{
				cntBadBeta++;
				continue;
			}
			
			if (OATHConst.isNA(tokens[KeyIdx[metaIdx][OATHConst.se]]))
			{
				cntBadSe++;
				continue;
			}

//			if (KeyIdx[metaIdx][FConstant.BP] != -1 && ConstValues.isNA(tokens[KeyIdx[metaIdx][FConstant.BP]]))
//			{
//				cntBadBp++;
//				continue;
//			}

			if (tokens[KeyIdx[metaIdx][OATHConst.refale]].length() != 1)
			{
				cntBadA1++;
				continue;
			}

			if (KeyIdx[metaIdx][OATHConst.altale] != -1 && tokens[KeyIdx[metaIdx][OATHConst.altale]].length() != 1)
			{
				cntBadA2++;
				continue;
			}

			if (ConstValues.isNA(tokens[KeyIdx[metaIdx][OATHConst.raf]]))
			{
				cntBadFreq++;
				continue;
			}

			if (ConstValues.isNA(tokens[KeyIdx[metaIdx][OATHConst.vg]]))
			{
				cntBadVg++;
				continue;
			}

			SynthFStat ms = new SynthFStat(tokens[KeyIdx[metaIdx][OATHConst.snp]], Float.parseFloat(tokens[KeyIdx[metaIdx][OATHConst.vg]]), Float.parseFloat(tokens[KeyIdx[metaIdx][OATHConst.beta]]), Float.parseFloat(tokens[KeyIdx[metaIdx][OATHConst.se]]), tokens[KeyIdx[metaIdx][OATHConst.refale]].charAt(0), tokens[KeyIdx[metaIdx][OATHConst.altale]].charAt(0));
			if (KeyIdx[metaIdx][OATHConst.chr] != -1)
			{
				if(tokens[KeyIdx[metaIdx][OATHConst.chr]].equalsIgnoreCase("X"))
				{
					ms.setChr(23);
				}
				else if(tokens[KeyIdx[metaIdx][OATHConst.chr]].equalsIgnoreCase("Y"))
				{
					ms.setChr(24);
				}
				else if(tokens[KeyIdx[metaIdx][OATHConst.chr]].equalsIgnoreCase("XY"))
				{
					ms.setChr(25);
				}
				else if(tokens[KeyIdx[metaIdx][OATHConst.chr]].equalsIgnoreCase("MT"))
				{
					ms.setChr(26);
				}
				else
				{
					int chr = -1;
					try
					{
						chr = Integer.parseInt(tokens[KeyIdx[metaIdx][OATHConst.chr]]);						
					}
					catch (NumberFormatException e)
					{
						Logger.printUserLog(e.toString() + " in line " + total + " in '" + workingMetaFile.get(metaIdx) + ".' is a bad value for chromosome. Skipped this marker.");
						continue;
					}
					ms.setChr(chr);
				}

				if (isChr)
				{
					if( ms.getChr() != chrKeep)
					{
						continue;
					}
				}
			}
			if (KeyIdx[metaIdx][OATHConst.bp] != -1)
			{
				long bp = -1;
				try
				{
					bp = Long.parseLong(tokens[KeyIdx[metaIdx][OATHConst.bp]]);	
				}
				catch (NumberFormatException e)
				{
					Logger.printUserLog(e.toString() + " in line " + total + " in '" + workingMetaFile.get(metaIdx) + "' is a bad value for position. Skipped this marker.");
					continue;
				}
				ms.setBP(bp);
			}

			if(sumstat.containsKey(ms.getSNP()))
			{
				Logger.printUserLog("Warning: Marker '" + ms.getSNP() +"' duplicated input, first instance used, others skipped.");
				cntDup++;
			}
			else
			{
				sumstat.put(ms.getSNP(), ms);
				snpArray.add(ms.getSNP());

				if ( MetaSNPTable.containsKey(ms.getSNP()))
				{
					ArrayList<Integer> snpCnt = MetaSNPTable.get(ms.getSNP());
					snpCnt.set(metaIdx, 1);
					Integer Int = snpCnt.get(snpCnt.size()-1);
					Int++;
					snpCnt.set(snpCnt.size()-1, Int);
				}
				else
				{
					ArrayList<Integer> snpCnt = NewIt.newArrayList();
					snpCnt.ensureCapacity(workingMetaFile.size()+1);
					for(int ii = 0; ii < workingMetaFile.size() + 1; ii++)
					{
						snpCnt.add(0);
					}
					snpCnt.set(metaIdx, 1);
					snpCnt.set(snpCnt.size()-1, 1);
					MetaSNPTable.put(ms.getSNP(), snpCnt);
				}
				cnt++;
			}
		}
		reader.close();

		if (cntMissSNP > 0)
		{
			String lc = cntMissSNP == 1 ? "locus" : "loci";
			Logger.printUserLog("Removed " + cntMissSNP + " " + lc + " due to bad marker name(s)");
		}

		if (cntBadChr > 0)
		{
			String lc = cntBadChr == 1 ? "locus" : "loci";
			Logger.printUserLog("Removed " + cntBadChr + " " + lc + " due to incorrect chromosome(s).");
		}

//		if (cntBadBp > 0)
//		{
//			String lc = cntBadBp == 1 ? "locus" : "loci";
//			Logger.printUserLog("Removed " + cntBadBp + " " + lc + " due to incorrect physical position(s).");				
//		}

		if (cntBadFreq > 0)
		{
			String lc = cntBadFreq == 1 ? "locus" : "loci";
			Logger.printUserLog("Removed " + cntBadFreq + " " + lc + " due to incorrect frequency.");
		}

		if (cntBadVg > 0)
		{
			String lc = cntBadVg == 1 ? "locus" : "loci";
			Logger.printUserLog("Removed " + cntBadVg + " " + lc + " due to incorrect genotypic variance.");
		}

		if (cntBadBeta > 0)
		{
			String lc = cntBadBeta == 1 ? "locus" : "loci";
			Logger.printUserLog("Removed " + cntBadBeta + " " + lc + " due to incorrect effect(s).");
		}

		if (cntBadSe > 0)
		{
			String lc = cntBadSe == 1 ? "locus" : "loci";
			Logger.printUserLog("Removed " + cntBadSe + " " + lc + " due to incorrect se.");
		}

		if (cntBadA1 > 0)
		{
			String lc = cntBadA1 == 1 ? "locus" : "loci";
			Logger.printUserLog("Removed " + cntBadA1 + " " + lc + " due to bad a1 allele(s).");
		}

		if (cntBadA2 > 0)
		{
			String lc = cntBadA2 == 1 ? "locus" : "loci";
			Logger.printUserLog("Removed " + cntBadA2 + " " + lc + " due to bad a2 allele(s).");
		}

		if (cntDup > 0)
		{
			String lc = cntDup == 1 ? "locus" : "loci";
			Logger.printUserLog("Removed " + cntDup + " duplicated " + lc);
		}

		MetaSNPArray.add(snpArray);
		if (cnt == 0)
		{
			Logger.printUserLog("Did not find any summary statistics from '" + workingMetaFile.get(metaIdx)+ "'.");
			System.exit(0);
		}
		else
		{
			Logger.printUserLog("Read " + cnt +" (of " + total + ") summary statistics from '" + workingMetaFile.get(metaIdx) + "'.");
		}

		return sumstat;
	}

	public ArrayList<String> getWorkingMetaFile()
	{
		return workingMetaFile;
	}

	private boolean isChr;
	private int chrKeep=0;
	private String[] field;
	private boolean isGZ;
	private int[][] KeyIdx; //snp, chr, bp, beta, se, p, a1, a2
	private ArrayList<String> workingMetaFile;
	private ArrayList<HashMap<String, SynthFStat>> MStat = NewIt.newArrayList();
	private ArrayList<ArrayList<String>> MetaSNPArray = NewIt.newArrayList();
	private HashMap<String, ArrayList<Integer>> MetaSNPTable = NewIt.newHashMap();

//	private int[] keepCohortIdx;

}
