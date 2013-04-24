package realcheck;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

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

public class RealCheckOne
{
	private GenotypeMatrix G1;

	private int[] markerIdx;

	private ArrayList<SNP> snpList;

	private ArrayList<PersonIndex> PersonTable1;

	private SampleFilter sf1;

	public RealCheckOne()
	{
		PLINKParser pp1 = null;
		if (CmdArgs.INSTANCE.getBinaryDataArgs(0).isSet())
		{
			pp1 = new PLINKBinaryParser(CmdArgs.INSTANCE.getBinaryDataArgs(0)
					.getBed(), CmdArgs.INSTANCE.getBinaryDataArgs(0)
					.getBim(), CmdArgs.INSTANCE.getBinaryDataArgs(0)
					.getFam());
		} else
		{
			Logger.printUserError("--bfile is not set.");
			System.exit(1);
		}
		pp1.Parse();

		sf1 = new SampleFilter(pp1.getPedigreeData(), pp1.getMapData());
		G1 = new GenotypeMatrix(sf1.getSample());
		PersonTable1 = sf1.getSample();
		snpList = pp1.getMapData().getMarkerList();

	}

	public void Check()
	{

		StringBuffer sb = new StringBuffer();
		sb.append(CmdArgs.INSTANCE.out);
		sb.append(".real");
		PrintStream ps = FileProcessor.CreatePrintStream(sb.toString());

		if (CmdArgs.INSTANCE.getRealCheckParameter().getSnps() != null)
		{
			Logger.printUserLog("A similarity matrix is generated with real-check SNPs.");
			getSelectedMarker();
		} else
		{
			getRandomMarker();
		}

		ps.print("file1.famid1 file1.id1 file1.famid2 file1.id2 score nmiss/loci\n");
		for (int i = 0; i < G1.getGRow(); i++)
		{
			for (int j = i; j < G1.getGRow(); j++)
			{
				double[] s = similarityScore(i, j);
				if (s[0] > CmdArgs.INSTANCE.getRealCheckParameter()
						.getThresholdLower()
						&& s[0] <= CmdArgs.INSTANCE.getRealCheckParameter()
								.getThresholdUpper())
				{
					PersonIndex ps1 = PersonTable1.get(i);
					PersonIndex ps2 = PersonTable1.get(j);
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

			int g1 = G1.getAdditiveScore(idx1, idx);
			int g2 = G1.getAdditiveScore(idx2, idx);
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

		markerIdx = new int[sf1.getMapFile().getMarkerList().size()];

		for (int i = 0; i < markerIdx.length; i++)
			markerIdx[i] = i;

		StringBuffer sb = new StringBuffer();
		sb.append(CmdArgs.INSTANCE.out);
		sb.append(".realsnp");

		PrintStream ps = FileProcessor.CreatePrintStream(sb.toString());
		for (int i = 0; i < markerIdx.length; i++)
		{
			int idx = markerIdx[i];
			SNP snp = snpList.get(idx);
			ps.print(snp.getChromosome() + " " + snp.getName() + " "
					+ snp.getDistance() + " " + snp.getPosition() + " "
					+ snp.getFirstAllele() + " " + snp.getSecAllele() + "\n");
		}
		ps.close();
	}

	public void getRandomMarker()
	{
		int mn = 0;
		int nMarker = sf1.getMapFile().getMarkerList().size();
		if (CmdArgs.INSTANCE.getRealCheckParameter().getMarkerNumberFlag())
		{
			if (CmdArgs.INSTANCE.getRealCheckParameter().getMarkerNumber() > nMarker)
			{
				Logger.printUserLog("Real-check marker number is reduced to "
						+ nMarker + ".");
				mn = nMarker;
			} else
			{
				mn = CmdArgs.INSTANCE.getRealCheckParameter()
						.getMarkerNumber();
			}
			markerIdx = new int[mn];
			RandomDataImpl rd = new RandomDataImpl();
			rd.reSeed(CmdArgs.INSTANCE.seed);

			markerIdx = rd.nextPermutation(markerIdx.length, mn);
			Arrays.sort(markerIdx);
		} else
		{
			markerIdx = new int[snpList.size()];
			for (int i = 0; i < markerIdx.length; i++)
			{
				markerIdx[i] = i;
			}
		}

		StringBuffer sb = new StringBuffer();
		sb.append(CmdArgs.INSTANCE.out);
		sb.append(".realsnp");

		PrintStream ps = FileProcessor.CreatePrintStream(sb.toString());
		for (int i = 0; i < markerIdx.length; i++)
		{
			int idx = markerIdx[i];
			SNP snp = snpList.get(idx);
			ps.print(snp.getChromosome() + " " + snp.getName() + " "
					+ snp.getDistance() + " " + snp.getPosition() + " "
					+ snp.getFirstAllele() + " " + snp.getSecAllele() + "\n");
		}
		ps.close();
	}
}
