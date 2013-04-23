package gear.profile;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import family.pedigree.file.SNP;
import family.plink.PLINKBinaryParser;
import family.plink.PLINKParser;
import family.popstat.GenotypeMatrix;
import family.qc.rowqc.SampleFilter;
import gear.Parameter;
import gear.profile.struct.QScore;
import gear.profile.struct.ScoreUnit;
import gear.util.FileProcessor;
import gear.util.Logger;
import gear.util.NewIt;
import gear.util.SNPMatch;

public class RiskScore
{

	private GenotypeMatrix G1;

	private HashMap<String, ScoreUnit> Score = NewIt.newHashMap();

	private ArrayList<SNP> snpList1;
	private HashMap<String, QScore> QS = NewIt.newHashMap();
	private double[][] q_score_range;
	private String[] QRName;
	private boolean isQ = false;

	ArrayList<Integer> scoreCoding = NewIt.newArrayList();

	private SampleFilter sf1;

	public RiskScore()
	{
		Logger.printUserLog("Generating risk profile for genotypes.");

		initial();

		PLINKParser pp1 = null;
		if (Parameter.INSTANCE.getFileParameter().isSet())
		{
			pp1 = new PLINKParser(Parameter.INSTANCE.getFileParameter()
					.getPedFile(), Parameter.INSTANCE.getFileParameter()
					.getMapFile());
		} else if (Parameter.INSTANCE.getBfileParameter(0).isSet())
		{
			pp1 = new PLINKBinaryParser(Parameter.INSTANCE.getBfileParameter(0)
					.getBedFile(), Parameter.INSTANCE.getBfileParameter(0)
					.getBimFile(), Parameter.INSTANCE.getBfileParameter(0)
					.getFamFile());
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

	private void initial()
	{
		// read score file
		gear.util.BufferedReader scoreReader = new gear.util.BufferedReader(
				Parameter.INSTANCE.getProfileParameter().getScoreFile(), "score");
		ScoreUnit scoreUnit;
		while ((scoreUnit = ScoreUnit.getNextScoreUnit(scoreReader)) != null)
		{
			Score.put(scoreUnit.getSNP(), scoreUnit);
		}
		scoreReader.close();

		Logger.printUserLog("Number of predictors: " + Score.size());

		// read q score file and q range file
		String qScoreFile = Parameter.INSTANCE.getProfileParameter().getQScoreFile(),
			   qScoreRangeFile = Parameter.INSTANCE.getProfileParameter().getQScoreRangeFile();
		if (qScoreFile != null && qScoreRangeFile != null)
		{
			// q score file
			gear.util.BufferedReader qScoreReader = new gear.util.BufferedReader(qScoreFile, "q-score");
			QScore qScore;
			while ((qScore = QScore.getNextQScore(qScoreReader)) != null)
			{
				QS.put(qScore.getSNP(), qScore);
			}
			qScoreReader.close();

			if (QS.size() == 0)
			{
				Logger.printUserError("Nothing is selected in '" + qScoreFile + "'.");
				System.exit(1);
			} else
			{
				Logger.printUserLog("Number of q-scores: " + QS.size());
			}

			// q range file
			gear.util.BufferedReader qRangeReader = new gear.util.BufferedReader(qScoreRangeFile, "q-score-range");
			ArrayList<ArrayList<String>> QR = NewIt.newArrayList();
			while (true)
			{
				String[] tokens = qRangeReader.readTokens(3);
				if (tokens == null)
				{
					break;
				}
				ArrayList<String> qr = NewIt.newArrayList();
				qr.add(tokens[0]); // range label
				qr.add(tokens[1]); // lower bound
				qr.add(tokens[2]); // upper bound
				QR.add(qr);
			}
			qRangeReader.close();

			if (QR.isEmpty())
			{
				Logger.printUserError("Nothing is selected in '" + qScoreRangeFile + "'.");
			} else
			{
				Logger.printUserLog("Number of q-ranges: " + QR.size());
			}

			q_score_range = new double[QR.size()][2];
			QRName = new String[QR.size()];

			for (int i = 0; i < QR.size(); i++)
			{
				ArrayList<String> qr = QR.get(i);
				QRName[i] = qr.get(0);
				q_score_range[i][0] = Double.parseDouble(qr.get(1));
				q_score_range[i][1] = Double.parseDouble(qr.get(2));
			}

			isQ = true;
		}
	}

	public void makeProfile()
	{
		if (isQ)
		{
			multipleProfile();
		} else
		{
			singleProfile();
		}
	}

	public void multipleProfile()
	{
		int Total = 0;
		int monoLocus = 0;
		int ATGCLocus = 0;
		int[] CCSNP = new int[QRName.length];
		int[][] GCInd = new int[G1.getGRow()][QRName.length];
		int[] matchScheme = new int[5];

		double[][] riskProfile = new double[G1.getGRow()][QRName.length];
		// int[] CC = new int[G1.getGRow()];
		for (int i = 0; i < snpList1.size(); i++)
		{

			boolean[] qsL2Flag = new boolean[QRName.length];
			Arrays.fill(qsL2Flag, false);

			SNP snp = snpList1.get(i);
			char a1 = snp.getRefAllele();
			char a2 = snp.getSecAllele();
			boolean isATGC = SNPMatch.Confusion(a1, a2);

			ScoreUnit su = null;
			boolean isMatchRef = true;
			double sc = 0;
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
					if (!Parameter.INSTANCE.keepATGC())
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

				if (su.getRefAllele().compareTo(Character.toString(a1)) == 0)
				{
					isMatchRef = false;
					matchScheme[0]++;

				} else if (su.getRefAllele().compareTo(Character.toString(a2)) == 0)
				{
					isMatchRef = true;
					matchScheme[1]++;

				} else if (su.getRefAllele().compareTo(
						SNPMatch.Flip(Character.toString(a1))) == 0)
				{
					isMatchRef = false;
					matchScheme[2]++;

				} else if (su.getRefAllele().compareTo(
						SNPMatch.Flip(Character.toString(a2))) == 0)
				{
					isMatchRef = true;
					matchScheme[3]++;

				} else
				{
					isMatchRef = false;
					matchScheme[4]++;
					continue;
				}

				sc = su.getScore();
				if (Parameter.INSTANCE.getTranFunction() == gear.RegressionModel.LOGIT)
				{
					if (isMatchRef)
					{
						sc = Math.log(sc);
					} else
					{
						sc = -1 * Math.log(sc);
					}
				} else
				{
					if (!isMatchRef)
					{
						sc = -1 * sc;
					}
				}
			} else
			{// this snp is not in the predictor panel;
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
						riskProfile[j][k] += sc * G1.getAdditiveScore(j, i);
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
				} else
				{
					riskProfile[i][j] /= 2 * GCInd[i][j];
				}
			}
		}

		Logger.printUserLog("Number of SNPs mapped to the score file: " + Total);
		Logger.printUserLog("Number of monomorphic loci removed: " + monoLocus);
		Logger.printUserLog("Number of ATGC loci "
				+ (Parameter.INSTANCE.keepATGC() ? "detected: " : "removed: ")
				+ ATGCLocus);

		for (int i = 0; i < QRName.length; i++)
		{
			Logger.printUserLog("Number of SNPs mapped into " + QRName[i]
					+ ": " + CCSNP[i]);
		}

		for (int i = 0; i < 4; i++)
		{
			Logger.printUserLog("Number of SNPs matching Scheme " + (1 + i)
					+ ": " + matchScheme[i]);
		}

		StringBuffer sbim = new StringBuffer();
		sbim.append(Parameter.INSTANCE.out);
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
		int[] matchScheme = new int[5];

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
			char a1 = snp.getRefAllele();
			char a2 = snp.getSecAllele();
			boolean isATGC = SNPMatch.Confusion(a1, a2);

			ScoreUnit su = null;
			boolean isMatchRef = false;
			double sc = 0;
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
					if (!Parameter.INSTANCE.keepATGC())
					{
						continue;
					}
				}
				su = Score.get(snp.getName());
				if (su.isMissing())
				{
					continue;
				}

				if (su.getRefAllele().equals(Character.toString(a1)))
				{
					isMatchRef = false;
					matchScheme[0]++;
					ArrayList<String> s = s4.get(0);
					s.add(snp.getName());
				} else if (su.getRefAllele().equals(Character.toString(a2)))
				{
					isMatchRef = true;
					matchScheme[1]++;
					ArrayList<String> s = s4.get(1);
					s.add(snp.getName());
				} else if (su.getRefAllele().equals(
						SNPMatch.Flip(Character.toString(a1))))
				{
					isMatchRef = false;
					matchScheme[2]++;
					ArrayList<String> s = s4.get(2);
					s.add(snp.getName());

				} else if (su.getRefAllele().equals(
						SNPMatch.Flip(Character.toString(a2))))
				{
					isMatchRef = true;
					matchScheme[3]++;
					ArrayList<String> s = s4.get(3);
					s.add(snp.getName());

				} else
				{
					isMatchRef = false;
					matchScheme[4]++;
					continue;
				}
				CCSNP++;

				sc = su.getScore();
				if (Parameter.INSTANCE.getTranFunction() == gear.RegressionModel.LOGIT)
				{// logit s
					if (isMatchRef)
					{
						sc = Math.log(sc);
					} else
					{
						sc = -1 * Math.log(sc);
					}
				}
			} else
			{// this snp is not in the predictor panel;
				continue;
			}

			for (int j = 0; j < G1.getGRow(); j++)
			{
				if (G1.getAdditiveScore(j, i) != GenotypeMatrix.missing)
				{
					riskProfile[j] += sc * G1.getAdditiveScore(j, i);
					GCInd[j]++;
				}
			}
		}

		for (int i = 0; i < riskProfile.length; i++)
		{
			if (GCInd[i] == 0)
			{
				riskProfile[i] = 0;
			} else
			{
				riskProfile[i] /= 2 * GCInd[i];
			}
		}

		Logger.printUserLog("Number of SNPs mapped to the score file: " + Total);
		Logger.printUserLog("Number of monomorphic loci removed: " + monoLocus);
		Logger.printUserLog("Number of ATGC loci "
				+ (Parameter.INSTANCE.keepATGC() ? "detected: " : "removed: ")
				+ ATGCLocus);
		Logger.printUserLog("Number of SNP scores in the score file: " + CCSNP);

		for (int i = 0; i < 4; i++)
		{
			Logger.printUserLog("Number of SNPs matching Scheme " + (1 + i)
					+ ": " + matchScheme[i]);
		}

		StringBuffer sbim = new StringBuffer();
		sbim.append(Parameter.INSTANCE.out);
		sbim.append(".profile");
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
