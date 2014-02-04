package gear.sumstat;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;


import gear.CmdArgs;
import gear.ConstValues;
import gear.data.Person;
import gear.family.pedigree.PersonIndex;
import gear.family.pedigree.file.MapFile;
import gear.family.plink.PLINKBinaryParser;
import gear.family.plink.PLINKParser;
import gear.family.popstat.GenotypeMatrix;
import gear.family.qc.rowqc.SampleFilter;
import gear.sumstat.qc.rowqc.SumStatQC;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;

public class Inbreeding
{
	private GenotypeMatrix G;
	private int numMarker;
	private double[] maf;
	private MapFile snpMap;
	private SumStatQC ssQC;

	private double[] N;
	private String[][] indKeep;
	private HashMap<String, Integer> groupID = NewIt.newHashMap();
	private double[][] mafGroup;
	private int[] IndGroup;
	private double[][] w;
	private double[] Fst;

	private int missing = -1;

	// private ArrayList<String> GroupInfor = NewIt.newArrayList();

	public Inbreeding()
	{
		PLINKParser pp = null;
		if (CmdArgs.INSTANCE.getFileArgs().isSet())
		{
			pp = new PLINKParser(CmdArgs.INSTANCE.getFileArgs()
					.getPed(), CmdArgs.INSTANCE.getFileArgs()
					.getMap());
		} else if (CmdArgs.INSTANCE.getBFileArgs(0).isSet())
		{
			pp = new PLINKBinaryParser(CmdArgs.INSTANCE.getBFileArgs(0)
					.getBed(), CmdArgs.INSTANCE.getBFileArgs(0)
					.getBim(), CmdArgs.INSTANCE.getBFileArgs(0)
					.getFam());
		} else
		{
			Logger.printUserError("No input files.");
			System.exit(1);
		}
		pp.Parse();
		SampleFilter sf = new SampleFilter(pp.getPedigreeData(),
				pp.getMapData());
		ssQC = new SumStatQC(pp.getPedigreeData(), pp.getMapData(), sf);
		snpMap = pp.getMapData();
		GenotypeMatrix gm = new GenotypeMatrix(ssQC.getSample());
		setup(gm);
	}

	public void setup(GenotypeMatrix gm)
	{
		G = gm;
		numMarker = G.getNumMarker();
		readKeepFile();
		readGroupInfor();
	}

	public void CalculateFst()
	{
		//The definition for FST is as defined in Population Genetics page 98 by John Gillespie.
		//It ranges between 0~1.
		int[][] g = G.getG();
		maf = new double[numMarker];
		mafGroup = new double[numMarker][groupID.size()];
		w = new double[numMarker][groupID.size()];
		N = new double[numMarker];

		Fst = new double[numMarker];
		for (int i = 0; i < g.length; i++)
		{
			int idx = IndGroup[i];
			if (idx == missing)
				continue;
			for (int j = 0; j < numMarker; j++)
			{
				int[] c = G.getBiAlleleGenotype(i, j);
				if (c[0] == 1)
				{
					maf[j]++;
					mafGroup[j][idx]++;
				}
				if (c[1] == 1)
				{
					maf[j]++;
					mafGroup[j][idx]++;
				}
				if (c[1] != Person.MissingAlleleCode)
				{
					w[j][idx]++;
					N[j]++;
				}
			}
		}

		PrintWriter fstOut = null;
		try
		{
			fstOut = new PrintWriter(
					new String(CmdArgs.INSTANCE.out + ".fst"));
		} catch (IOException e)
		{
			e.printStackTrace();
		}
		fstOut.print("CHR\tSNP\tBP\tRefA" + "\t");
		for (int i = 0; i < groupID.size(); i++)
		{
			fstOut.print("prop" + (i + 1) + "\t" + "Freq" + (i + 1) + "\t"
					+ "NInd" + (i + 1) + "\t");
		}
		fstOut.println("Freq\tFst");

		for (int i = 0; i < numMarker; i++)
		{
			double f = 0;
			maf[i] = maf[i] / (2 * N[i]);
			for (int j = 0; j < groupID.size(); j++)
			{
				if (w[i][j] != 0)
				{
					mafGroup[i][j] = mafGroup[i][j] / (2 * w[i][j]);
				} else
				{
					mafGroup[i][j] = 0;
				}
				w[i][j] = w[i][j] / (N[i]);
			}

			fstOut.print(snpMap.getSNP(i).getChromosome() + "\t"
					+ snpMap.getMarkerName(i) + "\t"
					+ snpMap.getSNP(i).getPosition() + "\t"
					+ snpMap.getSNP(i).getFirstAllele() + "\t");
			for (int j = 0; j < groupID.size(); j++)
			{
				f += w[i][j] * 2 * mafGroup[i][j] * (1 - mafGroup[i][j]);
				fstOut.print(w[i][j] + "\t" + mafGroup[i][j] + "\t"
						+ (int) (w[i][j] * N[i]) + "\t");
			}
			fstOut.print(maf[i] + "\t");
			if (maf[i] != 0)
			{
				Fst[i] = 1 - f / (2 * maf[i] * (1 - maf[i]));
			}
			fstOut.println(Fst[i]);
		}
		fstOut.close();
	}

	public void readGroupInfor()
	{
		ArrayList<PersonIndex> pt = ssQC.getSample();
		IndGroup = new int[pt.size()];
		Arrays.fill(IndGroup, missing);
		int i = 0;
		for (Iterator<PersonIndex> e = pt.iterator(); e.hasNext();)
		{
			PersonIndex pi = e.next();
			for (int j = 0; j < indKeep[0].length; j++)
			{
				if (pi.getFamilyID().compareTo(indKeep[0][j]) == 0
						&& pi.getIndividualID().compareTo(indKeep[1][j]) == 0)
				{
					String g = indKeep[2][j];
					int idx = groupID.get(g).intValue();
					IndGroup[i] = idx;
				}
			}
			i++;
		}
	}

	private void readKeepFile()
	{
		BufferedReader reader = FileUtil
				.FileOpen(CmdArgs.INSTANCE.fst_file);
		String line = null;
		ArrayList<String> famList = NewIt.newArrayList();
		ArrayList<String> indList = NewIt.newArrayList();
		ArrayList<String> groupList = NewIt.newArrayList();
		int gID = 0;
		try
		{
			while ((line = reader.readLine()) != null)
			{
				line.trim();
				String[] l = line.split(ConstValues.WHITESPACE_DELIMITER);
				if (l.length < 3)
					continue;
				famList.add(l[0]);
				indList.add(l[1]);
				groupList.add(l[2]);
				if (!groupID.containsKey(l[2]))
				{
					groupID.put(l[2], gID);
					gID++;
				}
			}
		} catch (IOException e)
		{
			Logger.handleException(e,
					"An exception occurred when reading the FST file '"
							+ CmdArgs.INSTANCE.fst_file + "'.");
		}
		indKeep = new String[3][];
		indKeep[0] = (String[]) famList.toArray(new String[0]);
		indKeep[1] = (String[]) indList.toArray(new String[0]);
		indKeep[2] = (String[]) groupList.toArray(new String[0]);

		Set<String> k = groupID.keySet();
		String[] g = new String[k.size()];
		for (Iterator<String> e = k.iterator(); e.hasNext();)
		{
			String K = e.next();
			int idx = groupID.get(K).intValue();
			g[idx] = K;
		}

		PrintWriter fstGrp = null;
		try
		{
			fstGrp = new PrintWriter(new String(CmdArgs.INSTANCE.out
					+ ".fst.grp"));
		} catch (IOException e)
		{
			e.printStackTrace();
		}

		for (int i = 0; i < g.length; i++)
		{
			fstGrp.append(g[i] + " " + (i + 1) + "\n");
		}
		fstGrp.close();
	}
}
