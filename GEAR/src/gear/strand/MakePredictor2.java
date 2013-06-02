package gear.strand;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import family.pedigree.file.SNP;
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
import gear.util.stat.Z;
import gear.util.structure.Predictor2;

public class MakePredictor2
{
	private GenotypeMatrix G1;

	private int[][] comSNPIdx;
	private double[][] allelefreq1;
	private double[] N1;
	// private double[] N2;
	private ArrayList<Boolean> flag;

	private String[] title;
	private ArrayList<SNP> snpList1;
	private ArrayList<Predictor2> predictorList = NewIt.newArrayList();
	ArrayList<Integer> scoreCoding = NewIt.newArrayList();

	private SampleFilter sf1;

	public MakePredictor2()
	{
		readPredictor();

		PLINKParser pp1 = null;
		if (CmdArgs.INSTANCE.getBFileArgs(0).isSet())
		{
			pp1 = new PLINKBinaryParser(CmdArgs.INSTANCE.getBFileArgs(0)
					.getBed(), CmdArgs.INSTANCE.getBFileArgs(0)
					.getBim(), CmdArgs.INSTANCE.getBFileArgs(0)
					.getFam());
		} else
		{
			Logger.printUserError("--bfile is not set");
			System.exit(1);
		}
		pp1.Parse();

		sf1 = new SampleFilter(pp1.getPedigreeData(), pp1.getMapData());
		G1 = new GenotypeMatrix(sf1.getSample());
		snpList1 = sf1.getMapFile().getMarkerList();
	}

	public void BuildPredictor()
	{
		DecimalFormat fmt = new DecimalFormat("#.###E0");

		StringBuffer sb = new StringBuffer();
		sb.append(CmdArgs.INSTANCE.out);
		sb.append(".mergesnp");
		PrintStream ps = FileProcessor.CreatePrintStream(sb.toString());
		ps.append("SNP\tChr\tPos\tA1_1st\tA2_1st\tA1_2nd\tA2_2nd\tMAF_A1_1st\tMAF_A1_2nd\tFlip\tMerged\tP\tScheme\n");

		allelefreq1 = new double[G1.getNumMarker()][3];
		N1 = new double[G1.getNumMarker()];
		flag = NewIt.newArrayList();

		CalculateAlleleFrequency(G1, allelefreq1, N1);

		getCommonSNP(snpList1);

		int qualified_snp = 0;
		for (int i = 0; i < comSNPIdx[0].length; i++)
		{
			int scheme = 0;
			boolean ATGCLocus = false;
			boolean flip = false;

			SNP snp1 = snpList1.get(comSNPIdx[0][i]);
			Predictor2 maf2 = predictorList.get(comSNPIdx[1][i]);
			char a1_1 = snp1.getFirstAllele();
			char a1_2 = snp1.getSecAllele();
			char a2_1 = maf2.getA1();
			char a2_2 = maf2.getA2();

			double ref1 = allelefreq1[comSNPIdx[0][i]][0];
			double ref2 = maf2.getMAF();
			boolean f = true;

			if (SNPMatch.IsBiallelic(a1_1, a1_2, a2_1, a2_2))
			{
				if (a1_1 == a2_1)
				{// scheme1
					scheme = 1;
					if (SNPMatch.isAmbiguous(a1_1, a1_2))
					{
						ATGCLocus = true;
						if (ref1 < 0.5 && ref2 < 0.5)
						{
							if (ref1 < CmdArgs.INSTANCE.getMergeArgs()
									.getMafCutoff()
									&& ref2 < CmdArgs.INSTANCE
											.getMergeArgs().getMafCutoff())
							{
								f = true;
							} else
							{
								f = false;
							}
							flip = false;
							scoreCoding.add(0);

						} else if (ref1 < 0.5 && ref2 > 0.5)
						{
							if (ref1 < CmdArgs.INSTANCE.getMergeArgs()
									.getMafCutoff()
									&& ref2 > 1 - CmdArgs.INSTANCE
											.getMergeArgs().getMafCutoff())
							{
								f = true;
							} else
							{
								f = false;
							}
							ref2 = 1 - ref1;
							flip = true;
							scoreCoding.add(1);

						} else if (ref1 > 0.5 && ref2 < 0.5)
						{
							if (ref1 > 1 - CmdArgs.INSTANCE
									.getMergeArgs().getMafCutoff()
									&& ref2 < CmdArgs.INSTANCE
											.getMergeArgs().getMafCutoff())
							{
								f = true;
							} else
							{
								f = false;
							}
							ref2 = 1 - ref2;
							flip = true;
							scoreCoding.add(1);

						} else
						{
							if (ref1 > 1 - CmdArgs.INSTANCE
									.getMergeArgs().getMafCutoff()
									&& ref2 > 1 - CmdArgs.INSTANCE
											.getMergeArgs().getMafCutoff())
							{
								f = true;
							} else
							{
								f = false;
							}
							flip = false;
							scoreCoding.add(0);
						}
						// debug
						f = false;
					} else
					{
						flip = false;
						scoreCoding.add(0);
						f = true;
					}

				} else if (a1_1 == a2_2)
				{// scheme2
					scheme = 2;
					if (SNPMatch.isAmbiguous(a1_1, a1_2))
					{
						ATGCLocus = true;
						if (ref1 < 0.5 && (1 - ref2) < 0.5)
						{
							if (ref1 < CmdArgs.INSTANCE.getMergeArgs()
									.getMafCutoff()
									&& (1 - ref2) < CmdArgs.INSTANCE
											.getMergeArgs().getMafCutoff())
							{
								f = true;
							} else
							{
								f = false;
							}
							scoreCoding.add(1);
							ref2 = 1 - ref2;
							flip = true;
						} else if (ref1 < 0.5 && (1 - ref2) > 0.5)
						{
							if (ref1 < CmdArgs.INSTANCE.getMergeArgs()
									.getMafCutoff()
									&& (1 - ref2) > 1 - CmdArgs.INSTANCE
											.getMergeArgs().getMafCutoff())
							{
								f = true;
							} else
							{
								f = false;
							}
							flip = false;
							scoreCoding.add(0);

						} else if (ref1 > 0.5 && (1 - ref2) < 0.5)
						{
							if (ref1 > 1 - CmdArgs.INSTANCE
									.getMergeArgs().getMafCutoff()
									&& (1 - ref2) < CmdArgs.INSTANCE
											.getMergeArgs().getMafCutoff())
							{
								f = true;
							} else
							{
								f = false;
							}
							flip = false;
							scoreCoding.add(0);

						} else
						{
							if (ref1 > 1 - CmdArgs.INSTANCE
									.getMergeArgs().getMafCutoff()
									&& (1 - ref2) > 1 - CmdArgs.INSTANCE
											.getMergeArgs().getMafCutoff())
							{
								f = true;
							} else
							{
								f = false;
							}
							flip = true;
							scoreCoding.add(1);
							ref2 = 1 - ref2;

						}
						// debug
						f = false;
					} else
					{
						flip = true;
						ref2 = 1 - ref2;
						scoreCoding.add(1);
						f = false;
					}
				} else if (a1_1 == SNPMatch.Flip(a2_1))
				{// scheme3
					scheme = 3;
					flip = true;
					scoreCoding.add(0);
					f = true;
				} else if (a1_1 == SNPMatch.Flip(a2_2))
				{// scheme4
					scheme = 4;
					flip = true;
					ref2 = 1 - ref2;
					scoreCoding.add(1);
					f = true;
					// debug
					f = false;
				} else
				{// outlier
					scheme = 5;
					f = false;
					scoreCoding.add(0);
				}

				double p = Z.OddsRatioTestPvalueTwoTail(ref1, ref2,
						N1[comSNPIdx[0][i]], maf2.getNInd() * 2);
				if (p < CmdArgs.INSTANCE.getMergeArgs().getPCutoff())
				{
					f = false;
				}
				if (!CmdArgs.INSTANCE.keepATGC() && ATGCLocus)
				{
					f = false;
				}
				if (CmdArgs.INSTANCE.removeFlip() && flip)
				{
					f = false;
				}
				flag.add(f);
				if (f)
					qualified_snp++;

				ps.println(snpList1.get(comSNPIdx[0][i]).getName() + " "
						+ snpList1.get(comSNPIdx[0][i]).getChromosome() + " "
						+ snpList1.get(comSNPIdx[0][i]).getPosition() + " "
						+ a1_1 + " " + a1_2 + " " + a2_1 + " " + a2_2 + " "
						+ " " + fmt.format(ref1) + " "
						+ fmt.format(maf2.getMAF()) + " " + flip + " " + f
						+ " " + p + " scheme" + scheme);
			} else
			{
				flag.add(false);
				scoreCoding.add(0);
			}
		}

		ps.close();
		if (qualified_snp == 0)
		{
			Logger.printUserError("Common SNPs between the two SNP files: None");
			System.exit(1);
		} else
		{
			Logger.printUserLog("Common SNP(s) between the two SNP files: "
					+ qualified_snp);
		}

		Logger.printUserLog("flag " + flag.size() + ": snpCoding "
				+ scoreCoding.size());
		WritePredictor();
	}

	private void getCommonSNP(ArrayList<SNP> snplist1)
	{
		HashMap<String, Integer> SNPMap = NewIt.newHashMap();
		for (Iterator<SNP> e = snplist1.iterator(); e.hasNext();)
		{
			SNP snp = e.next();
			SNPMap.put(snp.getName(), 0);
		}

		int c = 0;
		HashMap<String, Integer> SNPMapList2 = NewIt.newHashMap();
		for (int i = 0; i < predictorList.size(); i++)
		{
			Predictor2 maf = predictorList.get(i);
			String snp_name = maf.getSNP();
			if (SNPMap.containsKey(snp_name))
			{
				SNPMap.put(snp_name, 1);
				SNPMapList2.put(snp_name, i);
				c++;
			} else
			{
				SNPMap.put(snp_name, 0);
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
		int idx1 = 0;
		for (int i = 0; i < snplist1.size(); i++)
		{
			SNP snp = snplist1.get(i);
			String snp_name = snp.getName();
			if (SNPMap.containsKey(snp_name)
					&& SNPMap.get(snp_name).intValue() == 1)
			{
				comSNPIdx[0][idx1] = i;
				comSNPIdx[1][idx1] = SNPMapList2.get(snp_name).intValue();
				idx1++;
			}
		}
		Logger.printUserLog("idx1 " + idx1);

	}

	public void CalculateAlleleFrequency(GenotypeMatrix G, double[][] frq,
			double[] n)
	{
		int[][] g = G.getG();
		for (int i = 0; i < g.length; i++)
		{
			for (int j = 0; j < G.getNumMarker(); j++)
			{
				int[] c = G.getBiAlleleGenotype(i, j);
				frq[j][c[0]]++;
				frq[j][c[1]]++;
			}
		}
		for (int i = 0; i < G.getNumMarker(); i++)
		{
			double w = frq[i][0] + frq[i][1];
			n[i] = frq[i][0] + frq[i][1] + frq[i][2];
			if (w > 0)
			{
				for (int j = 0; j < frq[i].length - 1; j++)
				{
					frq[i][j] /= w;
				}
				frq[i][2] /= n[i];
			} else
			{
				frq[i][2] = 1;
			}
		}
	}

	public void readPredictor()
	{

		BufferedReader reader = FileProcessor.FileOpen(CmdArgs.INSTANCE
				.getPredictorFile());
		String line;
		try
		{
			line = reader.readLine();
			int idx = 1;
			line = line.trim();
			title = line.split("\\s+");

			while ((line = reader.readLine()) != null)
			{
				line = line.trim();
				Predictor2 maf = new Predictor2(line, title.length, idx++);
				predictorList.add(maf);
			}
		} catch (IOException e)
		{
			Logger.handleException(e,
					"An exception occurred when parsing the predictor file.");
		}
	}

	public void WritePredictor()
	{
		StringBuffer sbim = new StringBuffer();
		sbim.append(CmdArgs.INSTANCE.out);
		sbim.append(".predictor");
		PrintStream predictorFile = FileProcessor.CreatePrintStream(sbim
				.toString());

		int NMiss = 0;
		for (int i = 0; i < comSNPIdx[0].length; i++)
		{
			if (!flag.get(i))
			{
				continue;
			}
			SNP snp = snpList1.get(comSNPIdx[0][i]);
			Predictor2 pd = predictorList.get(comSNPIdx[1][i]);
			if (ConstValues.isNA(pd.getField(CmdArgs.INSTANCE.getPredictorIdx())))
			{
				NMiss++;
				continue;
			} else
			{
				double s = Double.parseDouble(pd.getField(CmdArgs.INSTANCE
						.getPredictorIdx()));
				if (CmdArgs.INSTANCE.getTranFunction() == gear.RegressionModel.LINEAR)
				{
					if (scoreCoding.get(i).intValue() == 1)
					{
						s *= -1;
					}
				} else
				{
					if (s < 0)
					{
						s = -9;
					} else
					{
						if (scoreCoding.get(i).intValue() == 1)
						{
							s = Math.log(1 / s);
						} else
						{
							s = Math.log(s);
						}
					}
				}
				predictorFile.append(snp.getName() + "\t" + snp.getFirstAllele()
						+ "\t" + s + "\n");
			}
		}
		predictorFile.close();
		Logger.printUserLog("Write preditor to " + sbim.toString());
		Logger.printUserLog(NMiss
				+ " SNP(s) have missing values and were not printed.");
	}

}
