package gear.strand;

import java.io.DataOutputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import family.pedigree.PersonIndex;
import family.pedigree.file.SNP;
import family.pedigree.genotype.BPerson;
import family.plink.PLINKBinaryParser;
import family.plink.PLINKParser;
import family.popstat.GenotypeMatrix;
import family.qc.rowqc.SampleFilter;
import gear.CmdArgs;
import gear.util.BufferedReader;
import gear.util.FileProcessor;
import gear.util.Logger;
import gear.util.NewIt;
import gear.util.SNPMatch;
import gear.util.stat.Z;
import gear.util.structure.MAF;

public class Strand
{
	private GenotypeMatrix genoMat;

	private int[][] comSNPIdx;
	private double[][] allelefreq1;
	private double[] N1;
	private double[] N2;
	private ArrayList<Boolean> flag;

	private ArrayList<SNP> snpList;
	private ArrayList<MAF> mafList = new ArrayList<MAF>();

	private ArrayList<PersonIndex> personIndexes;
	ArrayList<Integer> snpCoding = new ArrayList<Integer>();

	private SampleFilter sampleFilter;

	private byte byte1 = 108;
	private byte byte2 = 27;
	private byte byte3 = 1;

	private DataOutputStream os = null;

	public Strand()
	{
		readStrand();

		PLINKParser plinkParser = null;
		if (CmdArgs.INSTANCE.getBFileArgs(0).isSet())
		{
			plinkParser = new PLINKBinaryParser(CmdArgs.INSTANCE.getBFileArgs(0).getBed(),
			                                    CmdArgs.INSTANCE.getBFileArgs(0).getBim(),
			                                    CmdArgs.INSTANCE.getBFileArgs(0).getFam());
		}
		else
		{
			Logger.printUserError("--bfile is not set.");
			System.exit(1);
		}
		plinkParser.Parse();

		sampleFilter = new SampleFilter(plinkParser.getPedigreeData(), plinkParser.getMapData());
		genoMat = new GenotypeMatrix(sampleFilter.getSample());
		personIndexes = sampleFilter.getSample();
		snpList = sampleFilter.getMapFile().getMarkerList();
	}

	public void merge()
	{
		DecimalFormat fmt = new DecimalFormat("#.###E0");

		StringBuffer sb = new StringBuffer();
		sb.append(CmdArgs.INSTANCE.out);
		sb.append(".mergesnp");
		PrintStream ps = FileProcessor.CreatePrintStream(sb.toString());
		ps.append("SNP\tChr\tPos\tA1_1st\tA2_1st\tA1_2nd\tA2_2nd\tMAF_A1_1st\tMAF_A1_2nd\tFlip\tMerged\tP\tScheme\n");

		StringBuffer sb1 = new StringBuffer();
		sb1.append(CmdArgs.INSTANCE.out);
		sb1.append(".mergebadsnp");
		PrintStream ps1 = FileProcessor.CreatePrintStream(sb.toString());
		ps1.append("SNP\tChr\tPos\tA1_1st\tA2_1st\tA1_2nd\tA2_2nd\tMAF_A1_1st\tMAF_A1_2nd\tFlip\tMerged\tP\tScheme\n");

		allelefreq1 = new double[genoMat.getNumMarker()][3];
		N1 = new double[genoMat.getNumMarker()];
		flag = NewIt.newArrayList();

		calculateAlleleFrequency(genoMat, allelefreq1, N1);

		getCommonSNP(snpList);

		int qualified_snp = 0;
		for (int i = 0; i < comSNPIdx[0].length; i++)
		{
			int scheme = 0;
			boolean ATGCLocus = false;
			boolean flip = false;

			SNP snp1 = snpList.get(comSNPIdx[0][i]);
			MAF maf2 = mafList.get(comSNPIdx[1][i]);
			char a1_1 = snp1.getFirstAllele();
			char a1_2 = snp1.getSecAllele();
			char a2_1 = maf2.getAllele1();
			char a2_2 = maf2.getAllele2();

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
							}
							else
							{
								f = false;
							}
							flip = false;
							snpCoding.add(0);

						}
						else if (ref1 < 0.5 && ref2 > 0.5)
						{
							if (ref1 < CmdArgs.INSTANCE.getMergeArgs()
									.getMafCutoff()
									&& ref2 > 1 - CmdArgs.INSTANCE
											.getMergeArgs().getMafCutoff())
							{
								f = true;
							}
							else
							{
								f = false;
							}
							ref2 = 1 - ref1;
							flip = true;
							snpCoding.add(1);

						}
						else if (ref1 > 0.5 && ref2 < 0.5)
						{
							if (ref1 > 1 - CmdArgs.INSTANCE
									.getMergeArgs().getMafCutoff()
									&& ref2 < CmdArgs.INSTANCE
											.getMergeArgs().getMafCutoff())
							{
								f = true;
							}
							else
							{
								f = false;
							}
							ref2 = 1 - ref2;
							flip = true;
							snpCoding.add(1);

						}
						else
						{
							if (ref1 > 1 - CmdArgs.INSTANCE
									.getMergeArgs().getMafCutoff()
									&& ref2 > 1 - CmdArgs.INSTANCE
											.getMergeArgs().getMafCutoff())
							{
								f = true;
							}
							else
							{
								f = false;
							}
							flip = false;
							snpCoding.add(0);
						}
						// debug
						f = false;
					}
					else
					{
						flip = false;
						snpCoding.add(0);
						f = true;
					}

				}
				else if (a1_1 == a2_2)
				{// scheme2
					scheme = 2;
					if (SNPMatch.isAmbiguous(a1_1, a1_2))
					{
						ATGCLocus = true;
						if (ref1 < 0.5 && (1 - ref2) < 0.5)
						{
							if (ref1 < CmdArgs.INSTANCE.getMergeArgs()
									.getMafCutoff()
									&& 1 - ref2 < CmdArgs.INSTANCE
											.getMergeArgs().getMafCutoff())
							{
								f = true;
							} else
							{
								f = false;
							}
							snpCoding.add(1);
							ref2 = 1 - ref2;
							flip = true;
						}
						else if (ref1 < 0.5 && (1 - ref2) > 0.5)
						{
							if (ref1 < CmdArgs.INSTANCE.getMergeArgs()
									.getMafCutoff()
									&& 1 - ref2 > 1 - CmdArgs.INSTANCE
											.getMergeArgs().getMafCutoff())
							{
								f = true;
							} else
							{
								f = false;
							}
							flip = false;
							snpCoding.add(0);

						}
						else if (ref1 > 0.5 && (1 - ref2) < 0.5)
						{
							if (ref1 > 1 - CmdArgs.INSTANCE
									.getMergeArgs().getMafCutoff()
									&& 1 - ref2 < CmdArgs.INSTANCE
											.getMergeArgs().getMafCutoff())
							{
								f = true;
							} else
							{
								f = false;
							}
							flip = false;
							snpCoding.add(0);

						}
						else
						{
							if (ref1 > 1 - CmdArgs.INSTANCE
									.getMergeArgs().getMafCutoff()
									&& 1 - ref2 > 1 - CmdArgs.INSTANCE
											.getMergeArgs().getMafCutoff())
							{
								f = true;
							} else
							{
								f = false;
							}
							flip = true;
							snpCoding.add(1);
							ref2 = 1 - ref2;

						}
						// debug
						f = false;
					}
					else
					{
						flip = true;
						ref2 = 1 - ref2;
						snpCoding.add(1);
						f = false;
					}
				}
				else if (a1_1 == SNPMatch.Flip(a2_1))
				{// scheme3
					scheme = 3;
					flip = true;
					snpCoding.add(0);
					f = true;
				}
				else if (a1_1 == SNPMatch.Flip(a2_2))
				{// scheme4
					scheme = 4;
					flip = true;
					ref2 = 1 - ref2;
					snpCoding.add(1);
					f = true;
					// debug
					f = false;
				}
				else
				{// outlier
					scheme = 5;
					f = false;
					snpCoding.add(0);
				}

				double p = Z.OddsRatioTestPvalueTwoTail(ref1, ref2, N1[comSNPIdx[0][i]], N2[comSNPIdx[1][i]]);
				if (p < CmdArgs.INSTANCE.getMergeArgs().getMafCutoff())
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

				ps.println(snpList.get(comSNPIdx[0][i]).getName() + " "
						+ snpList.get(comSNPIdx[0][i]).getChromosome() + " "
						+ snpList.get(comSNPIdx[0][i]).getPosition() + " "
						+ a1_1 + " " + a1_2 + " " + a2_1 + " " + a2_2 + " "
						+ " " + fmt.format(ref1) + " "
						+ fmt.format(maf2.getMAF()) + " " + flip + " " + f
						+ " " + p + " scheme" + scheme);
			} else
			{
				flag.add(false);
				snpCoding.add(0);
				ps1.println(snpList.get(comSNPIdx[0][i]).getName() + " "
						+ snpList.get(comSNPIdx[0][i]).getChromosome() + " "
						+ snpList.get(comSNPIdx[0][i]).getPosition() + " "
						+ a1_1 + " " + a1_2 + " " + a2_1 + " " + a2_2 + " "
						+ " " + fmt.format(ref1) + " "
						+ fmt.format(maf2.getMAF()) + " " + flip + " " + f
						+ " " + 1 + " scheme" + scheme);
			}
		}

		ps.close();
		ps1.close();
		if (qualified_snp == 0)
		{
			Logger.printUserError("Number of common SNPs between the two SNP files: None");
			System.exit(1);
		} else
		{
			Logger.printUserLog("Number of common SNPs between the two SNP files: " + qualified_snp);
		}

		Logger.printUserLog("flag " + flag.size() + ": snpCoding " + snpCoding.size());
		writeFile();
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
		for (int i = 0; i < mafList.size(); i++)
		{
			MAF maf = mafList.get(i);
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
			Logger.printUserError("Number of common SNPs between the two SNP files: None");
			System.exit(1);
		} else
		{
			Logger.printUserLog("Number of common SNPs between the two SNP files: "
					+ c);
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

	public void calculateAlleleFrequency(GenotypeMatrix G, double[][] frq, double[] n)
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
			n[i] = frq[i][0] + frq[i][1];
			if (w > 0)
			{
				for (int j = 0; j < frq[i].length - 1; j++)
				{
					frq[i][j] /= w;
				}
				frq[i][2] /= frq[i][0] + frq[i][1] + frq[i][2];
			} else
			{
				frq[i][2] = 1;
			}
		}
	}

	public void readStrand()
	{
		BufferedReader reader = BufferedReader.openTextFile(CmdArgs.INSTANCE.getStrandFile(), "strand");
		reader.readNonEmptyLine();
		MAF maf;
		while ((maf = MAF.next(reader)) != null)
		{
			mafList.add(maf);
		}
	}

	public void writeFile()
	{
		StringBuffer sbim = new StringBuffer();
		sbim.append(CmdArgs.INSTANCE.out);
		sbim.append(".bim");
		PrintStream pbim = FileProcessor.CreatePrintStream(sbim.toString());

		for (int i = 0; i < comSNPIdx[0].length; i++)
		{
			if (flag.get(i))
			{
				SNP snp = snpList.get(comSNPIdx[0][i]);
				pbim.append(snp.getChromosome() + "\t" + snp.getName() + "\t" +
				            snp.getDistance() + "\t" + snp.getPosition() + "\t" +
				            snp.getFirstAllele() + "\t" + snp.getSecAllele() + "\n");
			}
		}
		pbim.close();

		StringBuffer sfam = new StringBuffer();
		sfam.append(CmdArgs.INSTANCE.out);
		sfam.append(".fam");
		PrintStream pfam = FileProcessor.CreatePrintStream(sfam.toString());
		for (Iterator<PersonIndex> e = personIndexes.iterator(); e.hasNext();)
		{
			PersonIndex per = e.next();
			BPerson bp = per.getPerson();
			pfam.append(bp.getFamilyID() + "\t" + bp.getPersonID() + "\t" +
			            bp.getDadID() + "\t" + bp.getMomID() + "\t" +
					    bp.getGender() + "\t" + bp.getAffectedStatus() + "\n");
		}

		pfam.close();

		StringBuffer sbed = new StringBuffer();
		sbed.append(CmdArgs.INSTANCE.out);
		sbed.append(".bed");
		try
		{
			os = new DataOutputStream(new FileOutputStream(sbed.toString()));
		}
		catch (FileNotFoundException e)
		{
			Logger.printUserError("Cannot create file '" + sbed.toString() + "'.");
			Logger.printUserError("Exception Message: " + e.getMessage());
			System.exit(1);
		}

		try
		{
			os.writeByte(byte1);
			os.writeByte(byte2);
			os.writeByte(byte3);

			for (int i = 0; i < comSNPIdx[0].length; i++)
			{
				if (!flag.get(i))
				{
					continue;
				}

				int snpIdx = comSNPIdx[0][i];
				byte gbyte = 0;
				int idx = 0;

				int posByte = snpIdx >> BPerson.shift;
				int posBite = (snpIdx & 0xf) << 1;

				for (int j = 0; j < personIndexes.size(); j++)
				{
					PersonIndex pi = personIndexes.get(j);
					BPerson bp = pi.getPerson();
					byte g = bp.getOriginalGenotypeScore(posByte, posBite);
					if (snpCoding.get(i).intValue() == 1)
					{
						switch (g)
						{
						case 0:
							g = 3;
							break;
						case 2:
							g = 2;
							break;
						case 3:
							g = 0;
							break;
						default:
							g = 1;
							break; // missing
						}
					}

					gbyte <<= 2 * idx;
					gbyte |= g;
					idx++;

					if (j != (personIndexes.size() - 1))
					{
						if (idx == 4)
						{
							os.writeByte(gbyte);
							gbyte = 0;
							idx = 0;
						}
					}
					else
					{
						os.writeByte(gbyte);
					}
				}
			}

			os.close();
		}
		catch (IOException e)
		{
			Logger.handleException(e, "An exception occurred when writing the bed file '" + sbed.toString() + "'.");
		}
	}
}
