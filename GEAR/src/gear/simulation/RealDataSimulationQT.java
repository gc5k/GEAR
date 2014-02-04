package gear.simulation;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Random;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math.random.RandomDataImpl;
import org.apache.commons.math.stat.StatUtils;

import gear.CmdArgs;
import gear.data.Person;
import gear.family.pedigree.PersonIndex;
import gear.family.pedigree.file.SNP;
import gear.family.plink.PLINKBinaryParser;
import gear.family.plink.PLINKParser;
import gear.family.qc.rowqc.SampleFilter;
import gear.simulation.gm.RealDataSimulationGenotypeMatrix;
import gear.simulation.qc.rowqc.*;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;

public class RealDataSimulationQT
{
	private String casualLociFile = null;
	private int[] casualLociIdx = null;
	private double[] b = null;
	private double[][] bv = null;
	private double[][] y = null;
	private Random rnd;
	private PLINKParser pp = null;
	private SampleFilter sf = null;
	RealDataSimulationGenotypeMatrix GM;
	private int sampleSize = 0;

	public RealDataSimulationQT()
	{
		if (CmdArgs.INSTANCE.getFileArgs().isSet())
		{
			pp = new PLINKParser(CmdArgs.INSTANCE.getFileArgs()
					.getPed(), CmdArgs.INSTANCE.getFileArgs()
					.getMap());
		} 
		else if (CmdArgs.INSTANCE.getBFileArgs(0).isSet())
		{
			pp = new PLINKBinaryParser(CmdArgs.INSTANCE.getBFileArgs(0)
					.getBed(), CmdArgs.INSTANCE.getBFileArgs(0)
					.getBim(), CmdArgs.INSTANCE.getBFileArgs(0)
					.getFam());
		}
		else
		{
			Logger.printUserError("No input files.");
			System.exit(1);
		}
		pp.Parse();
		sf = new SampleFilter(pp.getPedigreeData(), pp.getMapData());

	}

	public void GenerateSample()
	{
		rnd = new Random(CmdArgs.INSTANCE.simuSeed);
		if (CmdArgs.INSTANCE.simuCasualLoci != null)
		{
			getCasualLoci();
		} 
		else
		{
			getRandomCasualLoci(CmdArgs.INSTANCE.simuRndCasualLoci);
		}
		RealDataSimulationQC rdSimuQC = new RealDataSimulationQC(
				pp.getPedigreeData(), pp.getMapData(), sf);
		GM = new RealDataSimulationGenotypeMatrix(rdSimuQC);
		sampleSize = GM.getGRow();

		bv = new double[CmdArgs.INSTANCE.simuRep][sampleSize];
		b = generateAddEffects(casualLociIdx.length);
		y = new double[CmdArgs.INSTANCE.simuRep][];
		for (int i = 0; i < CmdArgs.INSTANCE.simuRep; i++)
		{
			bv[i] = calculateBV(b);
			double v = StatUtils.variance(bv[i]);
			double se = Math.sqrt(v / (CmdArgs.INSTANCE.simuHsq)
					* (1 - CmdArgs.INSTANCE.simuHsq));
			y[i] = generateY(bv[i], se);
			double[] t = new double[sampleSize];
			System.arraycopy(y[i], 0, t, 0, sampleSize);
			Arrays.sort(t);
		}

		StringBuffer sb1 = new StringBuffer(CmdArgs.INSTANCE.out);
		sb1.append(".eff");
		PrintStream ps1 = FileUtil.CreatePrintStream(sb1.toString());
		ArrayList<SNP> snpList = pp.getMapData().getMarkerList();
		for (int i = 0; i < casualLociIdx.length; i++)
		{
			SNP snp = snpList.get(casualLociIdx[i]);
			ps1.append(snp.getName() + " " + snp.getFirstAllele() + " "
					+ snp.getSecAllele() + " " + b[i] + "\n");
		}
		ps1.close();

		StringBuffer sb2 = new StringBuffer(CmdArgs.INSTANCE.out);
		sb2.append(".phen");
		PrintStream ps2 = FileUtil.CreatePrintStream(sb2.toString());
		ArrayList<PersonIndex> personTable = rdSimuQC.getSample();

		for (int i = 0; i < sampleSize; i++)
		{
			PersonIndex pi = personTable.get(i);
			ps2.append(pi.getFamilyID() + " " + pi.getIndividualID() + " ");
			for (int j = 0; j < CmdArgs.INSTANCE.simuRep; j++)
			{
				ps2.append(y[j][i] + " ");
			}
			ps2.append("\n");
		}
		ps2.close();
	}

	private double[] generateY(double[] bv, double se)
	{
		double[] phe = new double[bv.length];
		for (int i = 0; i < phe.length; i++)
		{
			phe[i] += bv[i] + rnd.nextGaussian() * se;
		}
		return phe;
	}

	private double[] calculateBV(double[] ae)
	{
		double[] bv = new double[sampleSize];
		for (int i = 0; i < bv.length; i++)
		{
			for (int j = 0; j < casualLociIdx.length; j++)
			{
				int idx = casualLociIdx[j];
				int g = GM.getGenotypeScore(i, idx);

				if (g == Person.MissingGenotypeCode)
					continue; // leave it alone if it is missing

				bv[i] += ae[j] * GM.getGenotypeScore(i, idx);
			}
		}
		return bv;
	}

	private double[] generateAddEffects(int len)
	{
		double[] effect = new double[len];
		for (int i = 0; i < effect.length; i++)
		{
			effect[i] = rnd.nextGaussian();
		}
		return effect;
	}

	private void getRandomCasualLoci(int num)
	{
		if (num > 0)
		{
			ArrayList<SNP> snpList = pp.getMapData().getMarkerList();

			RandomDataImpl rd = new RandomDataImpl();
			rd.reSeed(CmdArgs.INSTANCE.simuSeed);
			casualLociIdx = rd.nextPermutation(snpList.size(), num);
		}
		else 
		{
			casualLociIdx = new int[pp.getMapData().getMarkerNumber()];
			for (int i = 0; i < casualLociIdx.length; i++) 
			{
				casualLociIdx[i] = i;
			}
		}
	}

	private void getCasualLoci()
	{
		ArrayList<String> cl = NewIt.newArrayList();
		ArrayList<SNP> snpList = pp.getMapData().getMarkerList();

		if (CmdArgs.INSTANCE.simuCasualLoci != null)
		{
			casualLociFile = CmdArgs.INSTANCE.simuCasualLoci;
			BufferedReader reader = FileUtil.FileOpen(casualLociFile);
			String line = null;
			try
			{
				while ((line = reader.readLine()) != null)
				{
					String[] l = line.split("\\s+");
					if (l.length < 1)
						continue;
					cl.add(l[0]);
				}
			} 
			catch (IOException e)
			{
				Logger.handleException(e,
						"An exception occurred when reading the casual-loci file '"
								+ casualLociFile + "'.");
			}
		}

		if (cl.size() == 0)
		{
			casualLociIdx = new int[snpList.size()];
			for (int i = 0; i < casualLociIdx.length; i++)
				casualLociIdx[i] = i;
		} 
		else
		{
			ArrayList<Integer> Idx = NewIt.newArrayList();
			HashSet<String> SS = NewIt.newHashSet();
			for (int i = 0; i < cl.size(); i++)
			{
				SS.add(cl.get(i));
			}
			for (int i = 0; i < snpList.size(); i++)
			{
				SNP snp = snpList.get(i);
				String rs = snp.getName();
				if (SS.contains(rs))
				{
					Idx.add(i);
				}
			}
			Integer[] A = Idx.toArray(new Integer[0]);
			casualLociIdx = ArrayUtils.toPrimitive(A);
		}
	}

}
