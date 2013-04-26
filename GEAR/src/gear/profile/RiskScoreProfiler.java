package gear.profile;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import family.pedigree.file.SNP;
import family.plink.PLINKBinaryParser;
import family.plink.PLINKParser;
import family.popstat.GenotypeMatrix;
import family.qc.rowqc.SampleFilter;
import gear.CmdArgs;
import gear.profile.struct.QScore;
import gear.profile.struct.ScoreUnit;
import gear.util.FileProcessor;
import gear.util.Logger;
import gear.util.NewIt;
import gear.util.SNPMatch;

public class RiskScoreProfiler extends ProfilerBase
{
	private GenotypeMatrix G1;

	private ArrayList<SNP> snpList1;

	ArrayList<Integer> scoreCoding = NewIt.newArrayList();

	private SampleFilter sf1;

	public RiskScoreProfiler()
	{
		PLINKParser pp1 = null;
		if (CmdArgs.INSTANCE.getFileArgs().isSet())
		{
			pp1 = new PLINKParser(CmdArgs.INSTANCE.getFileArgs().getPed(),
					CmdArgs.INSTANCE.getFileArgs().getMap());
		} else if (CmdArgs.INSTANCE.getBFileArgs(0).isSet())
		{
			pp1 = new PLINKBinaryParser(CmdArgs.INSTANCE.getBFileArgs(0)
					.getBed(), CmdArgs.INSTANCE.getBFileArgs(0).getBim(),
					CmdArgs.INSTANCE.getBFileArgs(0).getFam());
		} else
		{
			Logger.printUserError("Neither --file nor --bfile is set.");
			System.exit(1);
		}
		pp1.Parse();

		sf1 = new SampleFilter(pp1.getPedigreeData(), pp1.getMapData());
		G1 = new GenotypeMatrix(sf1.getSample());
		snpList1 = sf1.getMapFile().getMarkerList();

	}

	public void multipleProfile()
	{
		int Total = 0;
		int monoLocus = 0;
		int ATGCLocus = 0;
		int[] CCSNP = new int[QRName.length];
		int[][] GCInd = new int[G1.getGRow()][QRName.length];
		int[] matchScheme;
		if (CmdArgs.INSTANCE.greedyFlag)
		{
			matchScheme = new int[5];
		} else
		{
			matchScheme = new int[3];
		}

		double[][] riskProfile = new double[G1.getGRow()][QRName.length];
		// int[] CC = new int[G1.getGRow()];
		for (int i = 0; i < snpList1.size(); i++)
		{

			boolean[] qsL2Flag = new boolean[QRName.length];
			Arrays.fill(qsL2Flag, false);

			SNP snp = snpList1.get(i);
			char a1 = snp.getFirstAllele();
			char a2 = snp.getSecAllele();
			boolean isATGC = SNPMatch.Confusion(a1, a2);

			ScoreUnit su = null;
			double sc = 0;
			int Tag = 0;
			if (Score.containsKey(snp.getName()))
			{
				Total++;
				if (!SNPMatch.IsBiallelic(a1, a2))
				{
					monoLocus++;
					continue;
				}

				if (isATGC)
				{
					ATGCLocus++;
					if (!CmdArgs.INSTANCE.keepATGC())
					{
						continue;
					}
				}
				su = Score.get(snp.getName());
				if (su.isMissing())
				{
					continue;
				}

				if (QS.containsKey(snp.getName()))
				{
					QScore qs = QS.get(snp.getName());
					if (!qs.isMissing())
					{
						for (int k = 0; k < QRName.length; k++)
						{
							if (qs.getQScore() >= q_score_range[k][0]
									&& qs.getQScore() <= q_score_range[k][1])
							{
								qsL2Flag[k] = true;
								CCSNP[k]++;
							}
						}
					}
				} else
				{
					continue;
				}

				if (!CmdArgs.INSTANCE.greedyFlag)
				{
					Tag = AsIs(su, a1, a2, matchScheme);
				} else
				{
					Tag = Greedy(su, a1, a2, matchScheme);
				}

				if (CmdArgs.INSTANCE.getTranFunction() == gear.RegressionModel.LOGIT)
				{// logit s
					sc = Math.log(su.getScore());
				} else
				{
					sc = su.getScore();
				}

			} else
			{// this snp is not in the predictor panel;
				continue;
			}

			if (Tag == ProfileConstant.MatchNeither) 
			{
				continue;
			}
			for (int k = 0; k < qsL2Flag.length; k++)
			{
				if (!qsL2Flag[k])
					continue;
				CCSNP[k]++;
				for (int j = 0; j < G1.getGRow(); j++)
				{
					if (G1.getAdditiveScore(j, i) != GenotypeMatrix.missing)
					{
						if (Tag == ProfileConstant.MatchRefAllele)
						{
							riskProfile[j][k] += sc * Model(G1.getAdditiveScore(j, i));
						} else if (Tag == ProfileConstant.MatchAltAllele)
						{
							riskProfile[j][k] += sc
								* Model((2 - G1.getAdditiveScore(j, i)));
						}
						GCInd[j][k]++;
					}
				}
			}
		}

		for (int i = 0; i < riskProfile.length; i++)
		{
			for (int j = 0; j < QRName.length; j++)
			{
				if (GCInd[i][j] == 0)
				{
					riskProfile[i][j] = 0;
				} else if (CmdArgs.INSTANCE.addFlag || CmdArgs.INSTANCE.domFlag || CmdArgs.INSTANCE.recFlag)
				{
					riskProfile[i][j] /= GCInd[i][j];
				} else
				{
					riskProfile[i][j] /= 2 * GCInd[i][j];
				}
			}
		}

		Logger.printUserLog("Number of SNPs mapped to the score file: " + Total);
		Logger.printUserLog("Number of monomorphic loci removed: " + monoLocus);
		Logger.printUserLog("Number of ATGC loci "
				+ (CmdArgs.INSTANCE.keepATGC() ? "detected: " : "removed: ")
				+ ATGCLocus);

		for (int i = 0; i < QRName.length; i++)
		{
			Logger.printUserLog("Number of SNPs mapped into " + QRName[i]
					+ ": " + CCSNP[i]);
		}

		for (int i = 0; i < matchScheme.length; i++)
		{
			Logger.printUserLog("Number of SNPs matching Scheme " + (1 + i)
					+ ": " + matchScheme[i]);
		}

		StringBuffer sbim = new StringBuffer();
		sbim.append(CmdArgs.INSTANCE.out);
		sbim.append(".profile");
		PrintStream predictorFile = FileProcessor.CreatePrintStream(sbim
				.toString());

		predictorFile.print("FID\tIID\tPHENO");
		for (int i = 0; i < QRName.length; i++)
		{
			predictorFile.print("\tScore." + QRName[i]);
		}
		predictorFile.println();
		for (int i = 0; i < riskProfile.length; i++)
		{
			predictorFile.print(sf1.getSample().get(i).getFamilyID() + "\t"
					+ sf1.getSample().get(i).getIndividualID() + "\t"
					+ sf1.getHukouBook().get(i).getCol6());
			for (int j = 0; j < riskProfile[i].length; j++)
			{
				predictorFile.print("\t" + riskProfile[i][j]);
			}
			predictorFile.println();
		}
		predictorFile.close();
	}

	public void singleProfile()
	{
		int Total = 0;
		int monoLocus = 0;
		int CCSNP = 0;
		int ATGCLocus = 0;
		int[] matchScheme;
		if (CmdArgs.INSTANCE.greedyFlag)
		{
			matchScheme = new int[5];
		} else
		{
			matchScheme = new int[3];
		}

		ArrayList<ArrayList<String>> s4 = NewIt.newArrayList();
		for (int i = 0; i < 4; i++)
		{
			ArrayList<String> s = NewIt.newArrayList();
			s4.add(s);
		}
		double[] riskProfile = new double[G1.getGRow()];
		int[] GCInd = new int[G1.getGRow()];
		for (int i = 0; i < snpList1.size(); i++)
		{
			SNP snp = snpList1.get(i);
			char a1 = snp.getFirstAllele();
			char a2 = snp.getSecAllele();
			boolean isATGC = SNPMatch.Confusion(a1, a2);

			ScoreUnit su = null;
			double sc = 0;
			int Tag = 0;

			if (Score.containsKey(snp.getName()))
			{
				Total++;
				if (!SNPMatch.IsBiallelic(a1, a2))
				{
					monoLocus++;
					continue;
				}

				if (isATGC)
				{
					ATGCLocus++;
					if (!CmdArgs.INSTANCE.keepATGC())
					{
						continue;
					}
				}
				su = Score.get(snp.getName());
				if (su.isMissing())
				{
					continue;
				}

				if (!CmdArgs.INSTANCE.greedyFlag)
				{
					Tag = AsIs(su, a1, a2, matchScheme);
				} else
				{
					Tag = Greedy(su, a1, a2, matchScheme);
				}
				CCSNP++;

				if (CmdArgs.INSTANCE.getTranFunction() == gear.RegressionModel.LOGIT)
				{// logit s
					sc = Math.log(su.getScore());
				} else
				{
					sc = su.getScore();
				}
			} else
			{// this snp is not in the predictor panel;
				continue;
			}

			if (Tag == ProfileConstant.MatchNeither)
			{
				continue;
			}
			
			for (int j = 0; j < G1.getGRow(); j++)
			{
				if (G1.getAdditiveScore(j, i) != GenotypeMatrix.missing)
				{
					if (Tag == ProfileConstant.MatchRefAllele)
					{
						riskProfile[j] += sc * Model(G1.getAdditiveScore(j, i));
					} else if (Tag == ProfileConstant.MatchAltAllele)
					{
						riskProfile[j] += sc
									* Model(2 - G1.getAdditiveScore(j, i));
					}
					GCInd[j]++;
				}
			}
		}

		for (int i = 0; i < riskProfile.length; i++)
		{
			if (GCInd[i] == 0)
			{
				riskProfile[i] = 0;
			} else if (CmdArgs.INSTANCE.addFlag || CmdArgs.INSTANCE.domFlag || CmdArgs.INSTANCE.recFlag)
			{
				riskProfile[i] /=  GCInd[i];
			} else 
			{
				riskProfile[i] /= 2 * GCInd[i];
			}
		}

		Logger.printUserLog("Number of SNPs mapped to the score file: " + Total);
		Logger.printUserLog("Number of monomorphic loci removed: " + monoLocus);
		Logger.printUserLog("Number of ATGC loci "
				+ (CmdArgs.INSTANCE.keepATGC() ? "detected: " : "removed: ")
				+ ATGCLocus);
		Logger.printUserLog("Number of SNP scores in the score file: " + CCSNP);

		for (int i = 0; i < matchScheme.length; i++)
		{
			Logger.printUserLog("Number of SNPs matching Scheme " + (1 + i)
					+ ": " + matchScheme[i]);
		}

		StringBuffer sbim = new StringBuffer();
		sbim.append(CmdArgs.INSTANCE.out);
		if (CmdArgs.INSTANCE.addFlag)
		{
			Logger.printUserLog("Additive model is used.");
			sbim.append(".add");
		} else if (CmdArgs.INSTANCE.domFlag) 
		{
			Logger.printUserLog("Dominant model is used.");
			sbim.append(".dom");
		} else if (CmdArgs.INSTANCE.recFlag)
		{
			Logger.printUserLog("Recessive model is used.");
			sbim.append(".rec");
		} else {
			Logger.printUserLog("Allelic model is used.");
		}
		sbim.append(".profile");
		Logger.printUserLog("Writing the profile into '" + sbim.toString() + "'.");
		PrintStream predictorFile = FileProcessor.CreatePrintStream(sbim
				.toString());
		predictorFile.println("FID\tIID\tPHENO\tSCORE");
		for (int i = 0; i < riskProfile.length; i++)
		{
			predictorFile.println(sf1.getSample().get(i).getFamilyID() + "\t"
					+ sf1.getSample().get(i).getIndividualID() + "\t"
					+ sf1.getHukouBook().get(i).getCol6() + "\t"
					+ riskProfile[i]);
		}
		predictorFile.close();
	}

}
