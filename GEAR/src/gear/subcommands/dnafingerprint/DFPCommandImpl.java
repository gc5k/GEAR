package gear.subcommands.dnafingerprint;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.TreeSet;

import org.apache.commons.math.random.RandomDataImpl;

import gear.data.Person;
import gear.family.pedigree.PersonIndex;
import gear.family.pedigree.file.SNP;
import gear.family.plink.PLINKBinaryParser;
import gear.family.plink.PLINKParser;
import gear.family.popstat.GenotypeMatrix;
import gear.family.qc.rowqc.SampleFilter;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;
import gear.util.SNPMatch;
import gear.util.pop.PopStat;

public class DFPCommandImpl extends CommandImpl
{
	private DFPCommandArguments dfpArgs;

	private GenotypeMatrix pGM1;
	private GenotypeMatrix pGM2;

	private double[][] allelefreq;
	private int[] markerIdx;

	private int[][] comSNPIdx;
	private HashMap<String, Integer> comSNPIdxMap;

	private ArrayList<SNP> refSNPList;

	private ArrayList<PersonIndex> PersonTable1;
	private ArrayList<PersonIndex> PersonTable2;

	private SampleFilter sf1;
	private SampleFilter sf2;

	private boolean[] snpMatch;
	
	@Override
	public void execute(CommandArguments cmdArgs)
	{
		dfpArgs = (DFPCommandArguments) cmdArgs;
		if (dfpArgs.getBFile2Flag() && dfpArgs.getBFile() != null)
		{
			RealCheck();
			Check();
		}
		else if (dfpArgs.getBFile() != null)
		{
			RealCheckOne();
			CheckOne();
		}
		else 
		{
			Logger.printUserError("--bfile or --bfile2 is not set.");
			System.exit(1);
		}
	}

	private void RealCheck()
	{
		PLINKParser pp1 = PLINKParser.parse((CommandArguments) dfpArgs);
		sf1 = new SampleFilter(pp1.getPedigreeData(), (CommandArguments) dfpArgs);
		pGM1 = new GenotypeMatrix(sf1.getSample(), pp1.getMapData(), (CommandArguments) dfpArgs);
		PersonTable1 = sf1.getSample();

		Logger.printUserLog("");
		Logger.printUserLog("Reading bfile2...");
		PLINKBinaryParser pp2 = new PLINKBinaryParser(dfpArgs.getBed2(), dfpArgs.getBim2(), dfpArgs.getFam2());
		pp2.Parse();
		sf2 = new SampleFilter(pp2.getPedigreeData());
		pGM2 = new GenotypeMatrix(sf2.getSample(), pp2.getMapData());
		PersonTable2 = sf2.getSample();
		Logger.printUserLog("");

		refSNPList = pGM1.getSNPList();
	}

	public void Check()
	{
		allelefreq = PopStat.calAlleleFrequency(pGM1);

		StringBuffer sb = new StringBuffer();
		sb.append(dfpArgs.getOutRoot());
		sb.append(".real");
		PrintStream ps = FileUtil.CreatePrintStream(sb.toString());

		getCommonSNP(pGM1.getSNPList(), pGM2.getSNPList());

		getRandomMarker();

		setSNPFlipFlag();

		double es = 0;
		double ss = 0;
		int n = 0;
		int identical = 0;

		ps.print("FID1 ID1 FID2 ID2 Match ExpMatch Score nmiss\n");
		for (int i = 0; i < pGM1.getGRow(); i++)
		{
			for (int j = 0; j < pGM2.getGRow(); j++)
			{
				double[] s = similarityScore(i, j);
				double ES = 0;
				double OS = 0;
				if (s[1] > 0)
				{
					ES = s[2]/s[1];
					OS = (s[0] - s[2]) / (s[1] - s[2]);
				}
				if (OS >= dfpArgs.getLowCutoff()
						&& OS <= dfpArgs.getHighCutoff())
				{
					PersonIndex ps1 = PersonTable1.get(i);
					PersonIndex ps2 = PersonTable2.get(j);
					ps.print(ps1.getFamilyID() + " " + ps1.getIndividualID()
								+ " " + ps2.getFamilyID() + " "
								+ ps2.getIndividualID() + " " + s[0] + " " + (ES * s[1]) + " " + OS + " " + s[1]
								+ "\n");
					identical++;
				}

				es += OS;
				ss += OS * OS;
				n++;

			}
		}
		ps.close();

		double E = 0;
		double v = 0;
		if (n > 0 )
		{
			E = es/n;
			v = Math.sqrt(ss/n - E * E);
		}
		else
		{
			E = 0;
			v = 0;
		}

		long N = pGM1.getGRow() * pGM2.getGRow();
		Logger.printUserLog("In total " + N + " individual pairs were compared.\n");
		Logger.printUserLog("Mean is: " + E);
		Logger.printUserLog("Standard deviation is: " + v);
//		Logger.printUserLog("Mean and SD were calculated with the exclusion of the pair of the identical individual.\n");
		
		double[] sChart = similarityScoreChart();
		Logger.printUserLog("=====Reference similarity score chart=====");
		Logger.printUserLog("Parent-offspring: " + (sChart[0] - sChart[3])/(1-sChart[3]));
		Logger.printUserLog("Full sib: " + (sChart[1] - sChart[3])/(1-sChart[3]) + "\n");
//		Logger.printUserLog("Half sib: " + sChart[2] + "\n");
		Logger.printUserLog(identical + " pairs were captured.");
		Logger.printUserLog("The result has been saved into '" + sb.toString() + "'.");
	}

	private double[] similarityScoreChart() {
		double[] sChart = new double[4];
		//0 for parent-offsprint
		//1 for full sib
		//2 for half sib
		//3 random
		for (int i = 0; i < markerIdx.length; i++)
		{
			int idx = markerIdx[i];
			double H = allelefreq[idx][1] * (1 - allelefreq[idx][1]);
			sChart[0] += 1-2*H;
			sChart[1] += (1-H) * (1-H) + 0.5 * H * H;
			sChart[2] += (1-H) * (1-H) + 0.5 * H * H;
			sChart[3] += 1 - 4 * H + 6 * H * H;
		}
		sChart[0] /= markerIdx.length;
		sChart[1] /= markerIdx.length;
		sChart[2] /= markerIdx.length;
		sChart[3] /= markerIdx.length;
		return sChart;
	}

	private double[] similarityScore(int idx1, int idx2)
	{
		double[] s = { 0, 0, 0 };

		for (int i = 0; i < markerIdx.length; i++)
		{

			int idx = markerIdx[i];

			int g1 = pGM1.getAdditiveScore(idx1, comSNPIdx[0][idx]);
			int g2 = pGM2.getAdditiveScore(idx2, comSNPIdx[1][idx]);
			if (g1 == Person.MissingGenotypeCode
					|| g2 == Person.MissingGenotypeCode)
				continue;
			if ( snpMatch[i] ) 
			{
				if (g1 == g2)
				{
					s[0]++;
				}
			}
			else 
			{
				if (g1 == (2-g2) )
				{
					s[0]++;
				}
			}
			s[1]++;
			double p = allelefreq[idx][1];
			s[2] += p * p * p * p + 4 * p * p * (1-p) * (1-p) + (1-p) * (1-p) * (1-p) * (1-p);

		}
		return s;
	}

	public void getRandomMarker()
	{
		int mn = 0;
		if (dfpArgs.getNumMarkerFlag())
		{
			if (dfpArgs.getNumMarker() > comSNPIdxMap
				.size())
			{
				Logger.printUserLog("Realcheck marker number was reduced to "
					+ comSNPIdxMap.size() + "\n");
				mn = comSNPIdxMap.size();
			} 
			else if (dfpArgs.getNumMarker() < 0)
			{
				mn = comSNPIdxMap.size();
			} else {
				mn = (int) dfpArgs.getNumMarker();
			}
		}
		else 
		{
			mn = comSNPIdxMap.size();
		}

		markerIdx = new int[mn];
		RandomDataImpl rd = new RandomDataImpl();
		rd.reSeed(dfpArgs.getSeed());

		markerIdx = rd.nextPermutation(comSNPIdxMap.size(), mn);

		Arrays.sort(markerIdx);
		StringBuffer sb = new StringBuffer();
		sb.append(dfpArgs.getOutRoot());
		sb.append(".realsnp");

		PrintStream ps = FileUtil.CreatePrintStream(sb.toString());
		for (int i = 0; i < markerIdx.length; i++)
		{
			int idx = markerIdx[i];
			SNP snp = refSNPList.get(comSNPIdx[0][idx]);
			ps.print(snp.getChromosome() + " " + snp.getName() + " "
					+ snp.getDistance() + " " + snp.getPosition() + " "
					+ snp.getFirstAllele() + " " + snp.getSecAllele() + "\n");
		}
		ps.close();
	}

	private void getCommonSNP(ArrayList<SNP> snplist1, ArrayList<SNP> snplist2)
	{
		HashMap<String, Boolean> SNPMap1 = NewIt.newHashMap();
		TreeSet<String> DupSNP1 = NewIt.newTreeSet();
		for (int i = 0; i < snplist1.size(); i++)
		{
			SNP snp = snplist1.get(i);
			if(SNPMap1.containsKey(snp.getName())) // removed duplicated snps
			{
				SNPMap1.remove(snp.getName());
				DupSNP1.add(snp.getName());
			}
			else
			{
				SNPMap1.put(snp.getName(), false);
			}
		}

		if (DupSNP1.size() > 0)
		{
			Logger.printUserLog(DupSNP1.size() + " SNPs with duplicated ID in bfile will not be used in real-check.");
		}
		if (SNPMap1.size() == 0)
		{
			Logger.printUserLog("No SNPs left for real-check.");
			Logger.printUserLog("GEAR quited.");
			System.exit(1);
		}

		HashMap<String, Integer> SNPMap2 = NewIt.newHashMap();
		TreeSet<String> DupSNP2 = NewIt.newTreeSet();
		for (int i = 0; i < snplist2.size(); i++)
		{
			SNP snp = snplist2.get(i);
			if(SNPMap2.containsKey(snp.getName())) //removed duplicated snps
			{
				SNPMap2.remove(snp.getName());
				DupSNP2.add(snp.getName());
			}
			else
			{
				SNPMap2.put(snp.getName(), i);
			}
		}

		if (DupSNP2.size() > 0)
		{
			Logger.printUserLog(DupSNP2.size() + " SNPs with duplicated ID in bfile2 will not be used in real-check.");
		}
		if (SNPMap2.size() == 0)
		{
			Logger.printUserLog("No SNPs left for real-check.");
			Logger.printUserLog("GEAR quited.");
			System.exit(1);
		}

		int c = 0;
		int ATGC = 0;
		HashMap<String, Integer> SNPMapList2 = NewIt.newHashMap();
		for (String key : SNPMap2.keySet())
		{
			String snp_name = key;
			int snpIdx = SNPMap2.get(key).intValue();
			SNP snp = snplist2.get(snpIdx);
			if (SNPMap1.containsKey(snp_name))
			{
				if (!SNPMatch.isAmbiguous(snp.getFirstAllele(), snp.getSecAllele())) 
				{
					SNPMap1.put(snp_name, true);
					SNPMapList2.put(snp_name, snpIdx);
					c++;
				}
				else 
				{
					ATGC++;
				}
			}
		}

		if (c == 0)
		{
			Logger.printUserError("Common SNPs between the two SNP files: None");
			System.exit(1);
		} 
		else
		{
			Logger.printUserLog("Common SNP(s) between the two SNP files: " + (c+ATGC));
			if (ATGC>0) 
			{
				Logger.printUserLog("Removed Ambiguous loci (A/T, G/C biallelic): " + ATGC);
			}
		}

		comSNPIdx = new int[2][c];
		comSNPIdxMap = NewIt.newHashMap();
		int idx1 = 0;
		for (int i = 0; i < snplist1.size(); i++)
		{
			SNP snp = snplist1.get(i);
			String snp_name = snp.getName();
			if (SNPMap1.containsKey(snp_name))
			{
				if (SNPMap1.get(snp_name).booleanValue())
				{
					comSNPIdxMap.put(snp.getName(), idx1);
					comSNPIdx[0][idx1] = i;
					comSNPIdx[1][idx1] = SNPMapList2.get(snp_name).intValue();
					idx1++;
				}
			}
		}
	}

	private void setSNPFlipFlag() 
	{
		snpMatch = new boolean[markerIdx.length];
		ArrayList<SNP> l1 = pGM1.getSNPList();
		ArrayList<SNP> l2 = pGM2.getSNPList();

		for (int i = 0; i < markerIdx.length; i++)
		{
			int idx = markerIdx[i];

			SNP s1 = l1.get(comSNPIdx[0][idx]);
			SNP s2 = l2.get(comSNPIdx[1][idx]);
			if(s1.getSecAllele() == s2.getSecAllele()) 
			{
				snpMatch[i] = true;
			}
			else if (s1.getSecAllele() == SNPMatch.Flip(s2.getSecAllele())) 
			{
				snpMatch[i] = true;
			}
			else
			{
				snpMatch[i] = false;
			}
		}
	}


////////////////////////////////////////////////
	
	public void RealCheckOne()
	{
		PLINKParser pp1 = PLINKParser.parse((CommandArguments) dfpArgs);
		sf1 = new SampleFilter(pp1.getPedigreeData(), (CommandArguments) dfpArgs);
		pGM1 = new GenotypeMatrix(sf1.getSample(), pp1.getMapData(), (CommandArguments) dfpArgs);
		PersonTable1 = sf1.getSample();
		refSNPList = pGM1.getSNPList();
		Logger.printUserLog("");
	}

	public void CheckOne()
	{

		allelefreq = PopStat.calAlleleFrequency(pGM1);

		StringBuffer sb = new StringBuffer();
		sb.append(dfpArgs.getOutRoot());
		sb.append(".real");
		PrintStream ps = FileUtil.CreatePrintStream(sb.toString());

		getRandomMarkerOne();

		StringBuffer sb1 = new StringBuffer();
		sb1.append(dfpArgs.getOutRoot());
		sb1.append(".realsnp");
		Logger.printUserLog(markerIdx.length + " realcheck SNPs have been saved into '" + sb1.toString() + "'.");

		double es = 0;
		double ss = 0;
		int n = 0;
		int identical = 0;
		ps.print("FID1 ID1 FID2 ID2 Match ExpMatch Score nmiss\n");
		for (int i = 0; i < pGM1.getGRow(); i++)
		{
			for (int j = 0; j < i; j++)
			{
				double[] s = similarityScoreOne(i, j);
				double ES = 0;
				double OS = 0;
				if (s[1] > 0)
				{
					ES = s[2]/s[1];
					OS = (s[0] - s[2]) / (s[1] - s[2]);
				}
				if (OS >= dfpArgs.getLowCutoff()
						&& OS <= dfpArgs.getHighCutoff())
				{
					PersonIndex ps1 = PersonTable1.get(i);
					PersonIndex ps2 = PersonTable1.get(j);
					ps.print(ps1.getFamilyID() + " " + ps1.getIndividualID()
							+ " " + ps2.getFamilyID() + " "
							+ ps2.getIndividualID() + " " + s[0] + " " + (ES * s[1]) + " " + OS + " " + s[1]
							+ "\n");
					identical++;
				}
				if(i != j) 
				{
					es += OS;
					ss += OS * OS;
					n++;
				}
			}
		}
		ps.close();

		double E = 0;
		double v = 0;
		if (n > 0 )
		{
			E = es/n;
			v = Math.sqrt(ss/n - E * E);
		}
		else 
		{
			E = 0;
			v = 0;
		}
		
		long N = pGM1.getGRow() * (pGM1.getGRow() + 1)/2;
		Logger.printUserLog("In total " + N + " individual pairs were compared.\n");
		Logger.printUserLog("Mean is: " + E);
		Logger.printUserLog("Standard deviation is: " + v);
		Logger.printUserLog("Mean and SD were calculated with the exclusion of the pair of the individual.\n");
		
		double[] sChart = similarityScoreChartOne();
		Logger.printUserLog("=====Reference similarity score chart=====");
		Logger.printUserLog("Parent-offspring: " + (sChart[0] - sChart[3])/(1-sChart[3]));
		Logger.printUserLog("Full sib: " + (sChart[1] - sChart[3])/(1-sChart[3]) + "\n");
//		Logger.printUserLog("Half sib: " + sChart[2] + "\n");		
		Logger.printUserLog(identical + " pair(s) were captured.");
		Logger.printUserLog("The result has been saved into '" + sb.toString() + "'.");

	}

	private double[] similarityScoreChartOne() {
		double[] sChart = new double[4];
		//0 for parent-offsprint
		//1 for full sib
		//2 for half sib
		//4 random
		for (int i = 0; i < markerIdx.length; i++)
		{
			int idx = markerIdx[i];
			double H = allelefreq[idx][1] * (1 - allelefreq[idx][1]);
			double p = allelefreq[idx][1];
			sChart[0] += 1-2*H;
			sChart[1] += (1-H) * (1-H) + 0.5 * H * H;
			sChart[2] += (1-H) * (1-H) + 0.5 * H * H;
			sChart[3] += p * p * p * p + 4 * p * p * (1-p) * (1-p) + (1-p) * (1-p) * (1-p) * (1-p);
		}
		sChart[0] /= markerIdx.length;
		sChart[1] /= markerIdx.length;
		sChart[2] /= markerIdx.length;
		sChart[3] /= markerIdx.length;
		return sChart;
	}

	private double[] similarityScoreOne(int idx1, int idx2)
	{
		double[] s = { 0, 0, 0};

		for (int i = 0; i < markerIdx.length; i++)
		{

			int idx = markerIdx[i];

			int g1 = pGM1.getAdditiveScore(idx1, idx);
			int g2 = pGM1.getAdditiveScore(idx2, idx);
			if (g1 == Person.MissingGenotypeCode
					|| g2 == Person.MissingGenotypeCode)
				continue;
			if (g1 == g2)
			{
				s[0]++;
			}
			s[1]++;
			double p = allelefreq[idx][1];
			s[2] += p * p * p * p + 4 * p * p * (1-p) * (1-p) + (1-p) * (1-p) * (1-p) * (1-p);
		}

		return s;
	}

	public void getSelectedMarkerOne()
	{

		markerIdx = new int[pGM1.getSNPList().size()];

		for (int i = 0; i < markerIdx.length; i++)
			markerIdx[i] = i;

		StringBuffer sb = new StringBuffer();
		sb.append(dfpArgs.getOutRoot());
		sb.append(".realsnp");

		PrintStream ps = FileUtil.CreatePrintStream(sb.toString());
		for (int i = 0; i < markerIdx.length; i++)
		{
			int idx = markerIdx[i];
			SNP snp = refSNPList.get(idx);
			ps.print(snp.getChromosome() + " " + snp.getName() + " "
					+ snp.getDistance() + " " + snp.getPosition() + " "
					+ snp.getFirstAllele() + " " + snp.getSecAllele() + "\n");
		}
		ps.close();
	}

	public void getRandomMarkerOne()
	{
		int mn = 0;
		int nMarker = pGM1.getSNPList().size();
		if (dfpArgs.getNumMarkerFlag())
		{
			if (dfpArgs.getNumMarker() > nMarker)
			{
				Logger.printUserLog("Real-check marker number is reduced to "
						+ nMarker + ".");
				mn = nMarker;
			} 
			else
			{
				mn = (int) dfpArgs.getNumMarker();
			}
			markerIdx = new int[mn];
			RandomDataImpl rd = new RandomDataImpl();
			rd.reSeed(dfpArgs.getSeed());

			markerIdx = rd.nextPermutation(markerIdx.length, mn);
			Arrays.sort(markerIdx);
		} 
		else
		{
			markerIdx = new int[refSNPList.size()];
			for (int i = 0; i < markerIdx.length; i++)
			{
				markerIdx[i] = i;
			}
		}

		StringBuffer sb = new StringBuffer();
		sb.append(dfpArgs.getOutRoot());
		sb.append(".realsnp");

		PrintStream ps = FileUtil.CreatePrintStream(sb.toString());
		for (int i = 0; i < markerIdx.length; i++)
		{
			int idx = markerIdx[i];
			SNP snp = refSNPList.get(idx);
			ps.print(snp.getChromosome() + " " + snp.getName() + " "
					+ snp.getDistance() + " " + snp.getPosition() + " "
					+ snp.getFirstAllele() + " " + snp.getSecAllele() + "\n");
		}
		ps.close();
	}
}
