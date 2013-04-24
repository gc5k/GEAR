package realcheck;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;

import org.apache.commons.math.random.RandomDataImpl;

import family.pedigree.PersonIndex;
import family.pedigree.file.SNP;
import family.pedigree.genotype.BPerson;
import family.plink.PLINKBinaryParser;
import family.plink.PLINKParser;
import family.popstat.GenotypeMatrix;
import family.qc.rowqc.SampleFilter;
import gear.CmdArgs;
import gear.util.FileProcessor;
import gear.util.Logger;
import gear.util.NewIt;
import gear.util.TextHelper;

public class RealCheck
{
	private GenotypeMatrix G1;
	private GenotypeMatrix G2;

	private int[] markerIdx;

	private int[][] comSNPIdx;
	private HashMap<String, Integer> comSNPIdxMap;

	private ArrayList<SNP> snpList;

	private ArrayList<PersonIndex> PersonTable1;
	private ArrayList<PersonIndex> PersonTable2;

	private SampleFilter sf1;
	private SampleFilter sf2;

	public RealCheck()
	{
		PLINKParser pp1 = null;
		PLINKParser pp2 = null;
		if (CmdArgs.INSTANCE.getBinaryDataArgs(0).isSet()
				&& CmdArgs.INSTANCE.getBinaryDataArgs(1).isSet())
		{
			pp1 = new PLINKBinaryParser(CmdArgs.INSTANCE.getBinaryDataArgs(0)
					.getBedFile(), CmdArgs.INSTANCE.getBinaryDataArgs(0)
					.getBimFile(), CmdArgs.INSTANCE.getBinaryDataArgs(0)
					.getFamFile());

			pp2 = new PLINKBinaryParser(CmdArgs.INSTANCE.getBinaryDataArgs(1)
					.getBedFile(), CmdArgs.INSTANCE.getBinaryDataArgs(1)
					.getBimFile(), CmdArgs.INSTANCE.getBinaryDataArgs(1)
					.getFamFile());
		} else
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
		} else
		{
			getRandomMarker();
		}

		ps.print("file1.famid file1.id file2.famid file2.id score nmiss/loci\n");
		for (int i = 0; i < G1.getGRow(); i++)
		{
			for (int j = 0; j < G2.getGRow(); j++)
			{
				double[] s = similarityScore(i, j);
				if (s[0] > CmdArgs.INSTANCE.getRealCheckParameter()
						.getThresholdLower()
						&& s[0] <= CmdArgs.INSTANCE.getRealCheckParameter()
								.getThresholdUpper())
				{
					PersonIndex ps1 = PersonTable1.get(i);
					PersonIndex ps2 = PersonTable2.get(j);
					ps.print(ps1.getFamilyID() + " " + ps1.getIndividualID()
							+ " " + ps2.getFamilyID() + " "
							+ ps2.getIndividualID() + " " + s[0] + " " + s[1]
							+ "/" + markerIdx.length + "\n");
				}
			}
		}
		ps.close();

	}

	private double[] similarityScore(int idx1, int idx2)
	{
		double[] s = { 0, 0 };

		for (int i = 0; i < markerIdx.length; i++)
		{

			int idx = markerIdx[i];

			int g1 = G1.getAdditiveScore(idx1, comSNPIdx[0][idx]);
			int g2 = G2.getAdditiveScore(idx2, comSNPIdx[1][idx]);
			if (g1 == BPerson.MissingGenotypeCode
					|| g2 == BPerson.MissingGenotypeCode)
				continue;
			if (g1 == g2)
			{
				s[0]++;
			}
			s[1]++;
		}

		if (s[1] > 0)
		{
			s[0] = s[0] / s[1];
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
		if (CmdArgs.INSTANCE.getRealCheckParameter().getMarkerNumber() > comSNPIdxMap
				.size())
		{
			Logger.printUserLog("Realcheck marker number was reduced to "
					+ comSNPIdxMap.size() + "\n");
			mn = comSNPIdxMap.size();
		} else
		{
			mn = CmdArgs.INSTANCE.getRealCheckParameter().getMarkerNumber();
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
		HashMap<String, Integer> SNPMap = NewIt.newHashMap();
		HashMap<String, String> SNPRef = NewIt.newHashMap();
		for (Iterator<SNP> e = snplist1.iterator(); e.hasNext();)
		{
			SNP snp = e.next();
			SNPMap.put(snp.getName(), 0);
			SNPRef.put(snp.getName(), Character.toString(snp.getFirstAllele()));
		}

		int c = 0;
		HashMap<String, Integer> SNPMapList2 = NewIt.newHashMap();
		for (int i = 0; i < snplist2.size(); i++)
		{
			SNP snp = snplist2.get(i);
			String snp_name = snp.getName();
			if (SNPMap.containsKey(snp_name))
			{
				if (SNPRef.get(snp_name).compareTo(
						Character.toString(snp.getFirstAllele())) == 0)
				{
					SNPMap.put(snp_name, 1);
					SNPMapList2.put(snp_name, i);
					c++;
				}
			}
		}

		for (Iterator<SNP> e = snplist2.iterator(); e.hasNext();)
		{
			SNP snp = e.next();
			if (SNPMap.containsKey(snp.getName()))
			{
				if (SNPRef.get(snp.getName()).compareTo(
						Character.toString(snp.getFirstAllele())) == 0)
				{
					SNPMap.put(snp.getName(), 1);
					c++;
				}
			}
		}

		if (c == 0)
		{
			Logger.printUserError("Common SNPs between the two SNP files: None");
			System.exit(1);
		} else
		{
			Logger.printUserLog("Common SNP(s) between the two SNP files: " + c);
		}

		comSNPIdx = new int[2][c];
		comSNPIdxMap = NewIt.newHashMap();
		int idx1 = 0;
		for (int i = 0; i < snplist1.size(); i++)
		{
			SNP snp = snplist1.get(i);
			String snp_name = snp.getName();
			if (SNPMap.containsKey(snp_name)
					&& SNPMap.get(snp_name).intValue() > 0)
			{
				comSNPIdxMap.put(snp.getName(), idx1);
				comSNPIdx[0][idx1] = i;
				comSNPIdx[1][idx1] = SNPMapList2.get(snp_name).intValue();
				idx1++;
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
				String[] l = line.split(TextHelper.WHITESPACE_DELIMITER);
				for (int i = 0; i < l.length; i++)
				{
					selectedSNP.add(l[i]);
				}
			}
		} catch (IOException e)
		{
			Logger.handleException(e,
					"An exception occurred when reading the real-check SNPs.");
		}
		Logger.printUserLog(selectedSNP.size() + " marker(s) is read in "
				+ CmdArgs.INSTANCE.getRealCheckParameter().getSnps() + ".");
		return selectedSNP;
	}
}
