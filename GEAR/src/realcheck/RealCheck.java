package realcheck;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import org.apache.commons.math.random.RandomDataImpl;

import family.pedigree.PersonIndex;
import family.pedigree.file.SNP;
import family.pedigree.genotype.BPerson;
import family.plink.PLINKBinaryParser;
import family.plink.PLINKParser;
import family.popstat.GenotypeMatrix;
import family.qc.rowqc.SampleFilter;
import gear.CmdArgs;
import gear.ConstValues;
import gear.util.FileProcessor;
import gear.util.Logger;
import gear.util.NewIt;
import gear.util.SNPMatch;
import gear.util.pop.PopStat;

public class RealCheck
{
	private GenotypeMatrix G1;
	private GenotypeMatrix G2;

	private double[][] allelefreq;
	private int[] markerIdx;

	private int[][] comSNPIdx;
	private HashMap<String, Integer> comSNPIdxMap;

	private ArrayList<SNP> snpList;

	private ArrayList<PersonIndex> PersonTable1;
	private ArrayList<PersonIndex> PersonTable2;

	private SampleFilter sf1;
	private SampleFilter sf2;

	private boolean[] snpMatch;
	public RealCheck()
	{
		PLINKParser pp1 = null;
		PLINKParser pp2 = null;
		if (CmdArgs.INSTANCE.getBFileArgs(0).isSet()
				&& CmdArgs.INSTANCE.getBFileArgs(1).isSet())
		{
			pp1 = new PLINKBinaryParser(CmdArgs.INSTANCE.getBFileArgs(0)
					.getBed(), CmdArgs.INSTANCE.getBFileArgs(0)
					.getBim(), CmdArgs.INSTANCE.getBFileArgs(0)
					.getFam());

			pp2 = new PLINKBinaryParser(CmdArgs.INSTANCE.getBFileArgs(1)
					.getBed(), CmdArgs.INSTANCE.getBFileArgs(1)
					.getBim(), CmdArgs.INSTANCE.getBFileArgs(1)
					.getFam());
		}
		else
		{
			Logger.printUserError("--bfile or --bfile2 is not set.");
			System.exit(1);
		}
		pp1.Parse();
		pp2.Parse();

		sf1 = new SampleFilter(pp1.getPedigreeData(), pp1.getMapData());
		sf2 = new SampleFilter(pp2.getPedigreeData(), pp2.getMapData());
		G1 = new GenotypeMatrix(sf1.getSample());
		G2 = new GenotypeMatrix(sf2.getSample());
		PersonTable1 = sf1.getSample();
		PersonTable2 = sf2.getSample();
		snpList = pp1.getMapData().getMarkerList();

	}

	public void Check()
	{
		allelefreq = PopStat.calAlleleFrequency(G1, snpList.size());

		StringBuffer sb = new StringBuffer();
		sb.append(CmdArgs.INSTANCE.out);
		sb.append(".real");
		PrintStream ps = FileProcessor.CreatePrintStream(sb.toString());

		getCommonSNP(sf1.getMapFile().getMarkerList(), sf2.getMapFile()
				.getMarkerList());

		if (CmdArgs.INSTANCE.getRealCheckParameter().getSnps() != null)
		{
			Logger.printUserLog("A similarity matrix is generated with real-check SNPs.");
			getSelectedMarker();
		}
		else
		{
			getRandomMarker();
		}

		setSNPFlipFlag();

		double es = 0;
		double ss = 0;
		int n = 0;

		ps.print("FID1 ID1 FID2 ID2 Match ExpMatch Score nmiss\n");
		for (int i = 0; i < G1.getGRow(); i++)
		{
			for (int j = 0; j < G2.getGRow(); j++)
			{
				double[] s = similarityScore(i, j);
				double ES = 0;
				double OS = 0;
				if (s[1] > 0)
				{
					ES = s[2]/s[1];
					OS = (s[0] - s[2]) / (s[1] - s[2]);
				}
				if (OS > CmdArgs.INSTANCE.getRealCheckParameter()
						.getThresholdLower()
						&& OS <= CmdArgs.INSTANCE.getRealCheckParameter()
								.getThresholdUpper())
				{
					PersonIndex ps1 = PersonTable1.get(i);
					PersonIndex ps2 = PersonTable2.get(j);
					ps.print(ps1.getFamilyID() + " " + ps1.getIndividualID()
								+ " " + ps2.getFamilyID() + " "
								+ ps2.getIndividualID() + " " + s[0] + " " + (ES * s[1]) + " " + OS + " " + s[1]
								+ "\n");
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

		long N = G1.getGRow() * G2.getGRow();
		Logger.printUserLog("In total " + N + " individual pairs were compared.\n");
		Logger.printUserLog("Mean is: " + E);
		Logger.printUserLog("Standard deviation is: " + v);
		Logger.printUserLog("Mean and SD were calculated with the exclusion of the pair of the individual.\n");
		
		double[] sChart = similarityScoreChart();
		Logger.printUserLog("=====Reference similarity score chart=====");
		Logger.printUserLog("Parent-offsprint: " + (sChart[0] - sChart[3])/(1-sChart[3]));
		Logger.printUserLog("Full sib: " + (sChart[1] - sChart[3])/(1-sChart[3]) + "\n");
//		Logger.printUserLog("Half sib: " + sChart[2] + "\n");		
		Logger.printUserLog("The result has been saved into '" + sb.toString() + "'.");
	}

	private double[] similarityScoreChart() {
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

	private double[] similarityScore(int idx1, int idx2)
	{
		double[] s = { 0, 0, 0 };

		for (int i = 0; i < markerIdx.length; i++)
		{

			int idx = markerIdx[i];

			int g1 = G1.getAdditiveScore(idx1, comSNPIdx[0][idx]);
			int g2 = G2.getAdditiveScore(idx2, comSNPIdx[1][idx]);
			if (g1 == BPerson.MissingGenotypeCode
					|| g2 == BPerson.MissingGenotypeCode)
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

	public void getSelectedMarker()
	{
		ArrayList<String> snps = readRealcheckSNPs();
		ArrayList<Integer> Idx = NewIt.newArrayList();
		for (int i = 0; i < snps.size(); i++)
		{
			String snp_name = snps.get(i);
			if (comSNPIdxMap.containsKey(snp_name))
			{
				Idx.add(comSNPIdxMap.get(snp_name));
			}
		}

		markerIdx = new int[Idx.size()];

		for (int i = 0; i < Idx.size(); i++)
			markerIdx[i] = Idx.get(i).intValue();
		Arrays.sort(markerIdx);

		StringBuffer sb = new StringBuffer();
		sb.append(CmdArgs.INSTANCE.out);
		sb.append(".realsnp");

		PrintStream ps = FileProcessor.CreatePrintStream(sb.toString());
		for (int i = 0; i < markerIdx.length; i++)
		{
			int idx = markerIdx[i];
			SNP snp = snpList.get(comSNPIdx[0][idx]);
			ps.print(snp.getChromosome() + " " + snp.getName() + " "
					+ snp.getDistance() + " " + snp.getPosition() + " "
					+ snp.getFirstAllele() + " " + snp.getSecAllele() + "\n");
		}
		ps.close();
	}

	public void getRandomMarker()
	{
		int mn = 0;
		if (CmdArgs.INSTANCE.getRealCheckParameter().getMarkerNumberFlag()) {
			if (CmdArgs.INSTANCE.getRealCheckParameter().getMarkerNumber() > comSNPIdxMap
				.size())
			{
				Logger.printUserLog("Realcheck marker number was reduced to "
					+ comSNPIdxMap.size() + "\n");
				mn = comSNPIdxMap.size();
			} 
			else
			{
				mn = CmdArgs.INSTANCE.getRealCheckParameter().getMarkerNumber();
			}
		}
		else 
		{
			mn = comSNPIdxMap.size();
		}

		markerIdx = new int[mn];
		RandomDataImpl rd = new RandomDataImpl();
		rd.reSeed(CmdArgs.INSTANCE.seed);

		markerIdx = rd.nextPermutation(comSNPIdxMap.size(), mn);

		Arrays.sort(markerIdx);
		StringBuffer sb = new StringBuffer();
		sb.append(CmdArgs.INSTANCE.out);
		sb.append(".realsnp");

		PrintStream ps = FileProcessor.CreatePrintStream(sb.toString());
		for (int i = 0; i < markerIdx.length; i++)
		{
			int idx = markerIdx[i];
			SNP snp = snpList.get(comSNPIdx[0][idx]);
			ps.print(snp.getChromosome() + " " + snp.getName() + " "
					+ snp.getDistance() + " " + snp.getPosition() + " "
					+ snp.getFirstAllele() + " " + snp.getSecAllele() + "\n");
		}
		ps.close();
	}

	private void getCommonSNP(ArrayList<SNP> snplist1, ArrayList<SNP> snplist2)
	{
		HashMap<String, Boolean> SNPMap = NewIt.newHashMap();
		for (int i = 0; i < snplist1.size(); i++)
		{
			SNP snp = snplist1.get(i);
			SNPMap.put(snp.getName(), false);
		}

		int c = 0;
		int ATGC = 0;
		HashMap<String, Integer> SNPMapList2 = NewIt.newHashMap();
		for (int i = 0; i < snplist2.size(); i++)
		{
			SNP snp = snplist2.get(i);
			String snp_name = snp.getName();
			if (SNPMap.containsKey(snp_name))
			{
				if (!SNPMatch.isAmbiguous(snp.getFirstAllele(), snp.getSecAllele())) 
				{
					SNPMap.put(snp_name, true);
					SNPMapList2.put(snp_name, i);
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
			Logger.printUserLog("Common SNP(s) between the two SNP files: " + c);
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
			if (SNPMap.get(snp_name).booleanValue())
			{
				comSNPIdxMap.put(snp.getName(), idx1);
				comSNPIdx[0][idx1] = i;
				comSNPIdx[1][idx1] = SNPMapList2.get(snp_name).intValue();
				idx1++;
			}
		}
	}

	private void setSNPFlipFlag() 
	{
		snpMatch = new boolean[markerIdx.length];
		ArrayList<SNP> l1 = sf1.getMapFile().getMarkerList();
		ArrayList<SNP> l2 = sf2.getMapFile().getMarkerList();

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

	private ArrayList<String> readRealcheckSNPs()
	{
		BufferedReader reader = FileProcessor.FileOpen(CmdArgs.INSTANCE
				.getRealCheckParameter().getSnps());
		String line = null;
		ArrayList<String> selectedSNP = NewIt.newArrayList();
		try
		{
			while ((line = reader.readLine()) != null)
			{
				String[] l = line.split(ConstValues.WHITESPACE_DELIMITER);
				for (int i = 0; i < l.length; i++)
				{
					selectedSNP.add(l[i]);
				}
			}
		} 
		catch (IOException e)
		{
			Logger.handleException(e,
					"An exception occurred when reading the real-check SNPs.");
		}
		Logger.printUserLog(selectedSNP.size() + " marker(s) is read in "
				+ CmdArgs.INSTANCE.getRealCheckParameter().getSnps() + ".");
		return selectedSNP;
	}
}
