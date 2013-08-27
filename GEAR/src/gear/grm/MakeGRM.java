package gear.grm;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;

import gear.HPC;
import gear.CmdArgs;
import gear.family.pedigree.Hukou;
import gear.family.pedigree.file.MapFile;
import gear.family.pedigree.file.PedigreeFile;
import gear.family.pedigree.file.SNP;
import gear.family.pedigree.genotype.BPerson;
import gear.family.plink.PLINKBinaryParser;
import gear.family.plink.PLINKParser;
import gear.family.popstat.GenotypeMatrix;
import gear.family.qc.rowqc.SampleFilter;
import gear.sumstat.qc.rowqc.SumStatQC;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;
import gear.util.pop.PopStat;
import gear.util.structure.MAF;

public class MakeGRM
{
	private double maf_threshold = 1e-5;

	private GenotypeMatrix G;
	private int numMarker;
	private double[][] allelefreq;
	private MapFile snpMap;
	private PedigreeFile pf;

	public MakeGRM()
	{
		PLINKParser pp = null;
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
		pf = pp.getPedigreeData();
		SampleFilter sf = new SampleFilter(pp.getPedigreeData(),
				pp.getMapData());
		SumStatQC ssQC = new SumStatQC(pp.getPedigreeData(), pp.getMapData(),
				sf);
		snpMap = pp.getMapData();
		GenotypeMatrix gm = new GenotypeMatrix(ssQC.getSample());
		G = gm;
		numMarker = G.getNumMarker();
		allelefreq = new double[numMarker][3];
	}

	public void makeGeneticRelationshipScore()
	{
		prepareMAF();
		StringBuffer sb = new StringBuffer();
		sb.append(CmdArgs.INSTANCE.out);
		PrintStream grm = null;
		BufferedWriter grmGZ = null;
		if (CmdArgs.INSTANCE.makeGRMTXTFlag)
		{
			sb.append(".grm.txt");
			grm = FileUtil.CreatePrintStream(sb.toString());
		} 
		else if (CmdArgs.INSTANCE.makeGRMFlag)
		{
			sb.append(".grm.gz");
			grmGZ = FileUtil.ZipFielWriter(sb.toString());
		}

		for (int i = 0; i < G.getGRow(); i++)
		{
			for (int j = 0; j <= i; j++)
			{
				double[] s = GRMScore(i, j);
				if (CmdArgs.INSTANCE.makeGRMTXTFlag)
				{
					grm.println((i + 1) + "\t" + (j + 1) + "\t" + s[0] + "\t"
							+ s[1]);
				} 
				else
				{
					try
					{
						grmGZ.append((i + 1) + "\t" + (j + 1) + "\t" + s[0]
								+ "\t" + s[1] + "\n");
					} 
					catch (IOException e)
					{
						Logger.handleException(e,
								"error in writing '" + sb.toString() + "' for "
										+ (i + 1) + " " + (j + 1) + ".");
					}
				}
			}
		}
		if (CmdArgs.INSTANCE.makeGRMTXTFlag)
		{
			grm.close();
		} 
		else
		{
			try
			{
				grmGZ.close();
			} 
			catch (IOException e)
			{
				Logger.handleException(e, " error in closing '" + sb.toString()
						+ "'.");
			}
		}
		Logger.printUserLog("Writing GRM scores into '" + sb.toString() + "'.");
		StringBuffer sb_id = new StringBuffer();
		sb_id.append(CmdArgs.INSTANCE.out);
		PrintStream grm_id = null;
		sb_id.append(".grm.id");
		grm_id = FileUtil.CreatePrintStream(sb_id.toString());

		ArrayList<Hukou> H = pf.getHukouBook();
		for (int i = 0; i < H.size(); i++)
		{
			Hukou h = H.get(i);
			grm_id.println(h.getFamilyID() + "\t" + h.getIndividualID());
		}
		grm_id.close();
		Logger.printUserLog("Writing individual information into '"
				+ sb_id.toString() + "'.");
	}

	public void GRMPartitioning(String[] args)
	{
		int N = G.getGRow();
		long T = N * (N + 1) / 2;

		int size = (int) Math.floor(T / CmdArgs.INSTANCE.grmPartition);

		for (int P = 1; P <= CmdArgs.INSTANCE.grmPartition; P++)
		{
			StringBuilder sb = new StringBuilder();
			sb.append(CmdArgs.INSTANCE.getHpcArgs().getName());
			sb.append("." + P + ".sh");

			PrintWriter pw = null;
			try
			{
				pw = new PrintWriter(sb.toString());
			} 
			catch (FileNotFoundException e)
			{
				Logger.handleException(e, "Cannot create the script file '"
						+ sb.toString() + "'.");
			}

			pw.println("#$ -cwd");
			int G = Integer.parseInt(CmdArgs.INSTANCE.getHpcArgs().getRam().substring(0, CmdArgs.INSTANCE.getHpcArgs().getRam().length()-1)) + 2;
			pw.println("#$ -l vf="
					+ G + "G");
			pw.println("#$ -l h_vmem="
					+ G + "G");

			StringBuilder nb = new StringBuilder();
			nb.append(CmdArgs.INSTANCE.getHpcArgs().getName() + "." + P);
			pw.println("#$ -N " + nb.toString());
			pw.println("#$ -m eas");
			pw.println("#$ -M "
					+ CmdArgs.INSTANCE.getHpcArgs().getEmail());

			pw.print("java -jar -Xmx"
					+ CmdArgs.INSTANCE.getHpcArgs().getRam() + " ");
			pw.print(HPC.class.getProtectionDomain().getCodeSource()
					.getLocation().getPath()
					+ " ");
			for (int i = 0; i < args.length; i++)
			{
				String arg = args[i];
				if (arg.equals("--shell") || arg.equals("--qsub"))
				{
					continue;
				}
				if (arg.equals("--email") || arg.equals("--ram")
						|| arg.equals("--name")
						|| arg.equals("--grm-partition")
						|| arg.equals("--grm_partition") || arg.equals("--out"))
				{
					++i;
					continue;
				}
				pw.print(arg + " ");
			}
			if (P != N)
			{
				pw.print("--grm-range " + (size * (P - 1) + 1) + ","
						+ (size * P) + " ");
			} 
			else
			{
				pw.print("--grm-range " + (size * (P - 1) + 1) + "," + T + " ");
			}
			pw.print("--out " + CmdArgs.INSTANCE.out + "." + P);
			pw.println();

			pw.close();
			Logger.printUserLog("Generated '" + sb.toString() + "'.");

			Runtime rt = Runtime.getRuntime();
			String cmd = "qsub " + sb.toString();
			Logger.printUserLog("Submitted '" + sb.toString() + "'.");
			try
			{
				rt.exec(cmd);
			} 
			catch (IOException e)
			{
				Logger.handleException(e, "Failed to execute command '" + cmd
						+ "'.");
			}
		}
	}

	public void makeGeneticRelationshipScore(int n0, int n1)
	{
		if (n0 < 1 || n1 > (G.getGRow() * G.getGRow()) / 2)
		{
			Logger.printUserError("incorrect range for grm pairs : " + n0 + ","
					+ n1);
			System.exit(0);
		} 
		else
		{
			Logger.printUserLog("generating grm scores for the pairs from "
					+ n0 + " to " + n1);
		}

		prepareMAF();
		StringBuffer sb = new StringBuffer();
		sb.append(CmdArgs.INSTANCE.out);
		PrintStream grm = null;
		BufferedWriter grmGZ = null;
		if (CmdArgs.INSTANCE.makeGRMTXTFlag)
		{
			sb.append(".grm.txt");
			grm = FileUtil.CreatePrintStream(sb.toString());
		} 
		else if (CmdArgs.INSTANCE.makeGRMFlag)
		{
			sb.append(".grm.gz");
			grmGZ = FileUtil.ZipFielWriter(sb.toString());
		}

		int i = 0, j = 0;
		for (; i < G.getGRow(); i++)
		{
			if ((i + 1) * (i + 1 + 1) / 2 >= n0)
				break;
		}
		j = n0 - i * (i + 1) / 2 - 1;
		int c = n0;

		while (c <= n1)
		{
			double[] s = GRMScore(i, j);
			if (CmdArgs.INSTANCE.makeGRMTXTFlag)
			{
				grm.println((i + 1) + "\t" + (j + 1) + "\t" + s[0] + "\t"
						+ s[1]);
			} 
			else
			{
				try
				{
					grmGZ.append((i + 1) + "\t" + (j + 1) + "\t" + s[0] + "\t"
							+ s[1] + "\n");
				} 
				catch (IOException e)
				{
					Logger.handleException(e,
							"error in writing '" + sb.toString() + "' for "
									+ (i + 1) + " " + (j + 1) + ".");
				}
			}
			if (j < i)
			{
				j++;
			} 
			else
			{
				i++;
				j = 0;
			}
			c++;
		}

		if (CmdArgs.INSTANCE.makeGRMTXTFlag)
		{
			grm.close();
		} 
		else
		{
			try
			{
				grmGZ.close();
			} 
			catch (IOException e)
			{
				Logger.handleException(e, " error in closing '" + sb.toString()
						+ "'.");
			}
		}
		Logger.printUserLog("Writing GRM scores into '" + sb.toString() + "'.");
		StringBuffer sb_id = new StringBuffer();
		sb_id.append(CmdArgs.INSTANCE.out);
		PrintStream grm_id = null;
		sb_id.append(".grm.id");
		grm_id = FileUtil.CreatePrintStream(sb_id.toString());

		ArrayList<Hukou> H = pf.getHukouBook();
		for (int k = 0; k < H.size(); k++)
		{
			Hukou h = H.get(k);
			grm_id.println(h.getFamilyID() + "\t" + h.getIndividualID());
		}
		grm_id.close();
		Logger.printUserLog("Writing individual information into '"
				+ sb_id.toString() + "'.");
	}

	private double[] GRMScore(int idx1, int idx2)
	{
		double[] s = { 0, 0 };

		for (int i = 0; i < allelefreq.length; i++)
		{

			if(allelefreq[i][1] == Double.NaN)
			{
				continue;
			}
			int g1 = G.getAdditiveScore(idx1, i);
			int g2 = G.getAdditiveScore(idx2, i);
			double m = allelefreq[i][1];
			if (g1 == BPerson.MissingGenotypeCode
					|| g2 == BPerson.MissingGenotypeCode)
			{
				continue;
			}
			else
			{
				if (m < maf_threshold || m > (1-maf_threshold))
				{// ignor too little maf
					continue;
				}
				else
				{
					s[0]++;
					s[1] += (g1 - 2 * m) * (g2 - 2 * m) / (2 * m * (1 - m));
				}
			}
		}

		if (s[0] > 0)
		{
			s[1] /= s[0];
		}
		else 
		{
			s[0] = 0;
			s[1] = 0;
		}

		return s;
	}

	private void prepareMAF()
	{
		if (CmdArgs.INSTANCE.ref_freq != null)
		{
			getRefFreq();

		} 
		else
		{
			allelefreq = PopStat.calAlleleFrequency(G, numMarker);
		}
	}

	private void getRefFreq()
	{
		HashMap<String, MAF> refMap = NewIt.newHashMap();
		gear.util.BufferedReader reader = gear.util.BufferedReader.openTextFile(CmdArgs.INSTANCE.ref_freq, "reference-frequency");
		reader.readNonEmptyLine();
		MAF maf;
		int numSnpsRead = 0;
		while ((maf = MAF.next(reader)) != null)
		{
			refMap.put(maf.getSNP(), maf);
			++numSnpsRead;
		}
		Logger.printUserLog("Read " + numSnpsRead + " SNPs in '" + CmdArgs.INSTANCE.ref_freq + "'.\n");

		ArrayList<SNP> snpList = snpMap.getMarkerList();
		int c = 0;
		for (int i = 0; i < snpList.size(); i++)
		{
			SNP snp = snpList.get(i);
			String snpName = snp.getName();
			double f = 0;
			if (refMap.containsKey(snpName))
			{
				f = refMap.get(snpName).getMAF();
			}
			if (allelefreq[i][1] <= 0.5)
			{
				if (allelefreq[i][1] < CmdArgs.INSTANCE.maf_range[0]
						|| allelefreq[i][1] > CmdArgs.INSTANCE.maf_range[1])
				{
					allelefreq[i][1] = 0;
				} 
				else
				{
					allelefreq[i][1] = f;
				}
			} 
			else
			{
				if ((1 - allelefreq[i][1]) < CmdArgs.INSTANCE.maf_range[0]
						|| (1 - allelefreq[i][1]) > CmdArgs.INSTANCE.maf_range[1])
				{
					allelefreq[i][1] = 0;
				} 
				else
				{
					allelefreq[i][1] = f;
				}
			}

			c++;
		}
		Logger.printUserLog("Got " + c + " matched reference alleles.\n");
	}

}
