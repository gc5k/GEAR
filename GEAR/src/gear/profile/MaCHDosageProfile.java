package gear.profile;

import gear.Parameter;
import gear.profile.struct.DosageInfor;
import gear.profile.struct.QScore;
import gear.profile.struct.ScoreUnit;
import gear.util.FileProcessor;
import gear.util.Logger;
import gear.util.NewIt;
import gear.util.SNPMatch;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

public class MaCHDosageProfile
{
	private String delim = "\\s+";
	private String[] dosageFile;
	private String[] inforFile;
	private HashMap<String, ScoreUnit> Score = NewIt.newHashMap();

	private ArrayList<String> ID = NewIt.newArrayList();

	private HashMap<String, QScore> QS = NewIt.newHashMap();
	private double[][] q_score_range;
	private String[] QRName;
	private boolean isQ = false;

	public MaCHDosageProfile()
	{
		Logger.printUserLog("Generating risk profiles for mach dosage.");
		initial();
	}

	private void initial()
	{
		// read score file
		String scoreFile = Parameter.INSTANCE.getProfileParameter().getScoreFile();
		if (scoreFile != null)
		{
			gear.util.BufferedReader scoreReader = new gear.util.BufferedReader(scoreFile, "score");
			ScoreUnit scoreUnit;
			while ((scoreUnit = ScoreUnit.getNextScoreUnit(scoreReader)) != null)
			{
				Score.put(scoreUnit.getSNP(), scoreUnit);
			}
		}

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

			if (QS.size() == 0)
			{
				Logger.printUserError("Nothing is selected in '" + qScoreFile + "'.");
				System.exit(1);
			} else
			{
				Logger.printUserLog("Number of SNP scores read from the q-score file '" + qScoreFile + "': " + QS.size());
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
				qr.add(tokens[0]); // label
				qr.add(tokens[1]); // lower bound
				qr.add(tokens[2]); // upper bound
				QR.add(qr);
			}

			if (QR.size() == 0)
			{
				Logger.printUserError("Nothing is selected in the q-range file '" + qScoreRangeFile + "'.");
			} else
			{
				Logger.printUserLog("Number of scores read from the q-range file '" + qScoreRangeFile + "': " + QR.size());
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

		String machDosageFile = Parameter.INSTANCE.getProfileParameter().getMachDosageFile();
		if (machDosageFile != null)
		{
			dosageFile = new String[1];
			dosageFile[0] = machDosageFile;
			File f = new File(dosageFile[0]);
			if (!f.exists())
			{
				Logger.printUserError("The dosage file '" + dosageFile[0] + "' does not exist");
				System.exit(1);
			}
			inforFile = new String[1];
			inforFile[0] = Parameter.INSTANCE.getProfileParameter().getMachInfoFile();
			f = new File(inforFile[0]);
			if (!f.exists())
			{
				Logger.printUserError("The information file " + inforFile[0] + "' does not exist");
				System.exit(1);
			}
		}
		else
		{
			String dosageBatchFile = Parameter.INSTANCE.getProfileParameter().getMachDosageBatchFile();
			if (dosageBatchFile == null)
			{
				Logger.printUserError("--mach-dosage-batch is not set.");
				System.exit(1);
			}
			String infoBatchFile = Parameter.INSTANCE.getProfileParameter().getMachInfoBatchFile();
			if (infoBatchFile == null)
			{
				Logger.printUserError("--mach-infor-batch is not set.");
				System.exit(1);
			}
			
			gear.util.BufferedReader dosageBatchReader = new gear.util.BufferedReader(dosageBatchFile, "MaCH dosage batch");
			ArrayList<String> l1 = NewIt.newArrayList();
			String line = null;
			while ((line = dosageBatchReader.readLine()) != null)
			{
				if (line.length() == 0)
					continue;
				l1.add(line);
			}
			
			dosageFile = (String[]) l1.toArray(new String[0]);
			for (int i = 0; i < dosageFile.length; i++)
			{
				File f = new File(dosageFile[i]);
				if (!f.exists())
				{
					Logger.printUserError("The dosage file '" + dosageFile[i] + "' does not exist.");
					System.exit(1);
				}
			}

			gear.util.BufferedReader infoBatchReader = new gear.util.BufferedReader(infoBatchFile, "MaCH information batch");
			ArrayList<String> l2 = NewIt.newArrayList();
			while ((line = infoBatchReader.readLine()) != null)
			{
				if (line.length() == 0)
					continue;
				l2.add(line);
			}
			
			inforFile = (String[]) l2.toArray(new String[0]);
			for (int i = 0; i < inforFile.length; i++)
			{
				File f = new File(inforFile[i]);
				if (!f.exists())
				{
					Logger.printUserError("The information file '" + inforFile[i] + "' does not exist.");
					System.exit(1);
				}
			}
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

		int monoLocus = 0;
		int[] CC = new int[QRName.length];
		int[] CCSNP = new int[QRName.length];
		double[][] riskProfile = null;
		int ATGCLocus = 0;

		int sumSNPMapped = 0;

		for (int i = 0; i < inforFile.length; i++)
		{

			ArrayList<DosageInfor> SD = readDosageInfor(inforFile[i]);

			ArrayList<ArrayList<Double>> dosage = readDosage(dosageFile[i], i);

			double[][] rs = new double[dosage.size()][QRName.length];

			if (i == 0)
				riskProfile = new double[dosage.size()][QRName.length];

			int[] c = new int[QRName.length];
			int cSNP = 0;

			for (int j = 0; j < SD.size(); j++)
			{

				boolean qsL1Flag = false;
				boolean[] qsL2Flag = new boolean[QRName.length];
				Arrays.fill(qsL2Flag, false);

				DosageInfor di = SD.get(j);
				String snp = di.getSNP();
				String refA = di.getRefAllele();
				String refB = di.getSecAllele();
				if (!SNPMatch.IsBiallelic(refA.charAt(0), refB.charAt(0)))
				{
					monoLocus++;
					continue;
				}

				if (di.isATGCLocus())
				{
					ATGCLocus++;
					if (!Parameter.INSTANCE.keepATGC())
					{
						continue;
					}
				}

				ScoreUnit su = null;
				if (Score.containsKey(snp))
				{
					su = Score.get(snp);
					if (su.isMissing())
					{
						continue;
					}
					cSNP++;

					if (QS.containsKey(snp))
					{
						QScore qs = QS.get(snp);
						if (!qs.isMissing())
						{
							qsL1Flag = true;
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
					}

				}
				else if (Parameter.INSTANCE.getProfileParameter().getScoreFile() == null)
				{
					su = new ScoreUnit(snp, refA, 1.0);
					cSNP++;
				}
				else
				{
					continue;
				}

				if (!qsL1Flag)
				{
					continue;
				}

				for (int k = 0; k < dosage.size(); k++)
				{

					double locusScore = 0;
					if (refA.compareTo(su.getRefAllele()) == 0)
					{
						locusScore += dosage.get(k).get(j).doubleValue()
								* su.getScore();
					} else if (refB.compareTo(su.getRefAllele()) == 0)
					{
						locusScore += (2 - dosage.get(k).get(j).doubleValue())
								* su.getScore();
					} else if (SNPMatch.Flip(refA).compareTo(su.getRefAllele()) == 0)
					{
						locusScore += dosage.get(k).get(j).doubleValue()
								* su.getScore();
					} else if (SNPMatch.Flip(refB).compareTo(su.getRefAllele()) == 0)
					{
						locusScore += (2 - dosage.get(k).get(j).doubleValue())
								* su.getScore();
					} else
					{
						continue;
					}

					for (int l = 0; l < qsL2Flag.length; l++)
					{
						if (!qsL2Flag[l])
							continue;
						rs[k][l] += locusScore;
						if (k == 0)
						{
							c[l]++;
							CC[l]++;
						}
					}
				}
			}

			sumSNPMapped += cSNP;
			Logger.printUserLog(dosageFile[i] + " mapped " + cSNP
					+ " SNP(s) to the score file.");
			for (int j = 0; j < c.length; j++)
			{
				Logger.printUserLog("\t" + c[j]
						+ " SNP(s) mapped to the range " + q_score_range[j][0]
						+ " " + q_score_range[j][1]);
			}

			for (int j = 0; j < riskProfile.length; j++)
			{
				for (int k = 0; k < riskProfile[j].length; k++)
				{
					riskProfile[j][k] += rs[j][k];
				}
			}
		}

		for (int i = 0; i < riskProfile.length; i++)
		{
			for (int j = 0; j < riskProfile[i].length; j++)
			{
				if (CC[j] == 0)
				{
					riskProfile[i][j] = 0;
				} else
				{
					riskProfile[i][j] /= 2 * CC[j];
				}
			}
		}

		for (int i = 0; i < CCSNP.length; i++)
		{
			Logger.printUserLog(CCSNP[i] + " SNP(s) mapped to the range "
					+ q_score_range[i][0] + " " + q_score_range[i][1]);
		}

		Logger.printUserLog("Number of monomorphic loci removed: " + monoLocus);
		Logger.printUserLog("Number of ATGC loci "
				+ (Parameter.INSTANCE.keepATGC() ? "detected: " : "removed: ")
				+ ATGCLocus);
		Logger.printUserLog("Number of SNPs mapped to the score file in total: "
				+ sumSNPMapped);

		StringBuffer sbim = new StringBuffer();
		sbim.append(Parameter.INSTANCE.out);
		sbim.append(".profile");
		PrintStream predictorFile = FileProcessor.CreatePrintStream(sbim
				.toString());
		predictorFile.print("FID\tIID");
		for (int i = 0; i < QRName.length; i++)
		{
			predictorFile.print("\tSCORE" + "." + QRName[i]);
		}
		predictorFile.println();

		for (int i = 0; i < riskProfile.length; i++)
		{
			String[] id = ID.get(i).split("->");
			predictorFile.print(id[0] + "\t" + id[1]);
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

		int monoLocus = 0;
		int CC = 0;
		int CCSNP = 0;
		double[] riskProfile = null;

		int ATGCLocus = 0;
		for (int i = 0; i < inforFile.length; i++)
		{

			ArrayList<DosageInfor> SD = readDosageInfor(inforFile[i]);

			ArrayList<ArrayList<Double>> dosage = readDosage(dosageFile[i], i);

			double[] rs = new double[dosage.size()];

			if (i == 0)
				riskProfile = new double[dosage.size()];

			int c = 0;
			int cSNP = 0;
			for (int j = 0; j < SD.size(); j++)
			{
				DosageInfor di = SD.get(j);
				String snp = di.getSNP();
				String refA = di.getRefAllele();
				String refB = di.getSecAllele();
				if (!SNPMatch.IsBiallelic(refA.charAt(0), refB.charAt(0)))
				{
					monoLocus++;
					continue;
				}

				if (di.isATGCLocus())
				{
					ATGCLocus++;
					if (Parameter.INSTANCE.keepATGC())
					{
						continue;
					}
				}

				ScoreUnit su = null;
				if (Score.containsKey(snp))
				{
					su = Score.get(snp);
					if (su.isMissing())
					{
						continue;
					}
					cSNP++;
					CCSNP++;
				} else
				{
					continue;
				}

				for (int k = 0; k < dosage.size(); k++)
				{
					if (refA.compareTo(su.getRefAllele()) == 0)
					{
						rs[k] += dosage.get(k).get(j).doubleValue()
								* su.getScore();
					} else if (refB.compareTo(su.getRefAllele()) == 0)
					{
						rs[k] += (2 - dosage.get(k).get(j).doubleValue())
								* su.getScore();
					} else if (SNPMatch.Flip(refA).compareTo(su.getRefAllele()) == 0)
					{
						rs[k] += dosage.get(k).get(j).doubleValue()
								* su.getScore();
					} else if (SNPMatch.Flip(refB).compareTo(su.getRefAllele()) == 0)
					{
						rs[k] += (2 - dosage.get(k).get(j).doubleValue())
								* su.getScore();
					} else
					{
						continue;
					}
					if (k == 0)
					{
						c++;
						CC++;
					}
				}
			}
			Logger.printUserLog(dosageFile[i] + " mapped " + cSNP
					+ " SNP(s) to the score file, and " + c
					+ " SNP(s) had scores.");

			for (int j = 0; j < rs.length; j++)
			{
				riskProfile[j] += rs[j];
			}
		}

		for (int i = 0; i < riskProfile.length; i++)
		{
			if (CC == 0)
			{
				riskProfile[i] = 0;
			} else
			{
				riskProfile[i] /= 2 * CC;
			}
		}

		Logger.printUserLog("Number of monomorphic loci removed: " + monoLocus);
		Logger.printUserLog("Number of ATGC loci "
				+ (Parameter.INSTANCE.keepATGC() ? "detected: " : "removed: ")
				+ ATGCLocus);
		Logger.printUserLog("Number of SNPs mapped to the score file in total: "
				+ CCSNP);
		Logger.printUserLog("Number of SNPs having scores: " + CC);

		StringBuffer sbim = new StringBuffer();
		sbim.append(Parameter.INSTANCE.out);
		sbim.append(".profile");
		PrintStream predictorFile = FileProcessor.CreatePrintStream(sbim
				.toString());
		predictorFile.println("FID\tIID\tSCORE");
		for (int i = 0; i < riskProfile.length; i++)
		{
			String[] id = ID.get(i).split("->");
			predictorFile.println(id[0] + "\t" + id[1] + "\t" + riskProfile[i]);
		}
		predictorFile.close();
	}

	private ArrayList<DosageInfor> readDosageInfor(String file)
	{
		Logger.printUserLog("Reading the dosage information file '" + file
				+ "'.");
		ArrayList<DosageInfor> SD = NewIt.newArrayList();

		BufferedReader readerDIFile = FileProcessor.FileOpen(file);
		String lineScore = null;

		try
		{
			lineScore = readerDIFile.readLine(); // header
			while ((lineScore = readerDIFile.readLine()) != null)
			{
				DosageInfor di = new DosageInfor(lineScore);
				SD.add(di);
			}
		} catch (IOException e)
		{
			Logger.handleException(e,
					"An exception occurred when reading the dosage information file '"
							+ file + "'.");
		}
		return SD;
	}

	private ArrayList<ArrayList<Double>> readDosage(String file, int idx)
	{
		Logger.printUserLog("Reading the dosage file '" + file + "'.");
		ArrayList<ArrayList<Double>> dosage = NewIt.newArrayList();
		BufferedReader readerDosageFile = FileProcessor.ZipFileOpen(file);
		String lineDosage = null;
		try
		{
			while ((lineDosage = readerDosageFile.readLine()) != null)
			{
				String s[] = lineDosage.split(delim);
				if (idx == 0)
				{
					ID.add(s[0]);
				}

				ArrayList<Double> d = NewIt.newArrayList();
				for (int i = 2; i < s.length; i++)
				{
					d.add(Double.parseDouble(s[i]));
				}
				dosage.add(d);
			}
		} catch (IOException e)
		{
			Logger.handleException(e,
					"An exception occurred when reading the dosage file '"
							+ file + "'.");
		}

		return dosage;
	}
}
