package gear.realcheck;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.math.random.RandomDataImpl;

import gear.CmdArgs;
import gear.family.pedigree.PersonIndex;
import gear.family.pedigree.file.SNP;
import gear.family.pedigree.genotype.BPerson;
import gear.family.plink.PLINKBinaryParser;
import gear.family.plink.PLINKParser;
import gear.family.popstat.GenotypeMatrix;
import gear.family.qc.rowqc.SampleFilter;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.pop.PopStat;

public class RealCheckOne
{
	private GenotypeMatrix G1;

	private double[][] allelefreq;
	private int[] markerIdx;

	private ArrayList<SNP> snpList;

	private ArrayList<PersonIndex> PersonTable1;

	private SampleFilter sf1;

	public RealCheckOne()
	{
		PLINKParser pp1 = null;
		if (CmdArgs.INSTANCE.getBFileArgs(0).isSet())
		{
			pp1 = new PLINKBinaryParser(CmdArgs.INSTANCE.getBFileArgs(0)
					.getBed(), CmdArgs.INSTANCE.getBFileArgs(0)
					.getBim(), CmdArgs.INSTANCE.getBFileArgs(0)
					.getFam());
		} 
		else
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

		allelefreq = PopStat.calAlleleFrequency(G1, snpList.size());

		StringBuffer sb = new StringBuffer();
		sb.append(CmdArgs.INSTANCE.out);
		sb.append(".real");
		PrintStream ps = FileUtil.CreatePrintStream(sb.toString());

		if (CmdArgs.INSTANCE.getRealCheckParameter().getSnps() != null)
		{
			getSelectedMarker();
		} 
		else
		{
			getRandomMarker();
		}

		StringBuffer sb1 = new StringBuffer();
		sb1.append(CmdArgs.INSTANCE.out);
		sb1.append(".realsnp");
		Logger.printUserLog(markerIdx.length + " realcheck SNPs have been saved into '" + sb1.toString() + "'.");

		double es = 0;
		double ss = 0;
		int n = 0;
		ps.print("FID1 ID1 FID2 ID2 Match ExpMatch Score nmiss\n");
		for (int i = 0; i < G1.getGRow(); i++)
		{
			for (int j = 0; j <= i; j++)
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
					PersonIndex ps2 = PersonTable1.get(j);
					ps.print(ps1.getFamilyID() + " " + ps1.getIndividualID()
							+ " " + ps2.getFamilyID() + " "
							+ ps2.getIndividualID() + " " + s[0] + " " + (ES * s[1]) + " " + OS + " " + s[1]
							+ "\n");
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
		
		long N = G1.getGRow() * (G1.getGRow() + 1)/2;
		Logger.printUserLog("In total " + N + " individual pairs were compared.\n");
		Logger.printUserLog("Mean is: " + E);
		Logger.printUserLog("Standard deviation is: " + v);
		Logger.printUserLog("Mean and SD were calculated with the exclusion of the pair of the individual.\n");
		
		double[] sChart = similarityScoreChart();
		Logger.printUserLog("=====Reference similarity score chart=====");
		Logger.printUserLog("Parent-offspring: " + (sChart[0] - sChart[3])/(1-sChart[3]));
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
		double[] s = { 0, 0, 0};

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
			double p = allelefreq[idx][1];
			s[2] += p * p * p * p + 4 * p * p * (1-p) * (1-p) + (1-p) * (1-p) * (1-p) * (1-p);
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

		PrintStream ps = FileUtil.CreatePrintStream(sb.toString());
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
			} 
			else
			{
				mn = CmdArgs.INSTANCE.getRealCheckParameter()
						.getMarkerNumber();
			}
			markerIdx = new int[mn];
			RandomDataImpl rd = new RandomDataImpl();
			rd.reSeed(CmdArgs.INSTANCE.simuSeed);

			markerIdx = rd.nextPermutation(markerIdx.length, mn);
			Arrays.sort(markerIdx);
		} 
		else
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

		PrintStream ps = FileUtil.CreatePrintStream(sb.toString());
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
