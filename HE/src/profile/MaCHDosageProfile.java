package profile;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import parameter.Parameter;
import profile.struct.DosageInfor;
import profile.struct.ScoreUnit;
import profile.struct.QScore;
import util.FileProcessor;
import util.NewIt;
import util.SNPMatch;

public class MaCHDosageProfile {
	private String delim = "\\s+";
	private Parameter par;
	private String[] dosageFile;
	private String[] inforFile;
	private String scoreFile;
	private HashMap<String, ScoreUnit> Score = NewIt.newHashMap();

	private ArrayList<String> ID = NewIt.newArrayList();

	private String q_score_file;
	private String q_score_range_file;

	private HashMap<String, QScore> QS = NewIt.newHashMap();
	private double[][] q_score_range;
	private String[] QRName;
	private boolean isQ = false;

	public MaCHDosageProfile (Parameter p) {
		par = p;
		initial();
	}

	private void initial() {

		//read score file
		scoreFile = par.scoreFile;
		BufferedReader readerScoreFile = FileProcessor.FileOpen(scoreFile);
		String lineScore = null;
		try {
			while((lineScore = readerScoreFile.readLine())!=null) {
				if (lineScore.length() == 0) continue;
				ScoreUnit su = new ScoreUnit(lineScore);
				Score.put(su.getSNP(), su);
			}
		} catch (IOException e) {
			e.printStackTrace(System.err);
			System.exit(0);
		}

		//read q score file and q range file
		if (par.q_score_file != null && par.q_score_range_file != null) {

			//q score file
			q_score_file = par.q_score_file;
			BufferedReader readerQScoreFile = FileProcessor.FileOpen(q_score_file);
			String lineQScore = null;
			try {
				while((lineQScore = readerQScoreFile.readLine()) != null) {
					if (lineQScore.length() == 0) continue;
					QScore qs = new QScore(lineQScore);
					QS.put(qs.getSNP(), qs);
				}
			} catch (IOException e) {
				e.printStackTrace(System.err);
				System.exit(0);
			}
			if (QS.size() == 0) {
				System.out.println("nothing has been selected in " + q_score_file);
				System.exit(0);
			} else {
				System.out.println("read in " + QS.size() + " SNP scores from " + q_score_file + ".");
			}

			//q range file
			q_score_range_file = par.q_score_range_file;
			BufferedReader readerQRangeFile = FileProcessor.FileOpen(q_score_range_file);
			String lineQRange = null;
			ArrayList<ArrayList<String>> QR = NewIt.newArrayList();
			try {
				while ((lineQRange = readerQRangeFile.readLine()) != null) {
					if (lineQRange.length() == 0) continue;
					String[] s = lineQRange.split(delim);
					if(s.length<3) continue;
					ArrayList<String> qr = NewIt.newArrayList();
					qr.add(s[0]); qr.add(s[1]); qr.add(s[2]);
					QR.add(qr);
				}
			} catch (IOException e) {
				e.printStackTrace(System.err);
				System.exit(0);
			}
			if (QR.size() == 0) {
				System.out.println("nothing has been selected in " + q_score_range_file);
			} else {
				System.out.println("read in " + QR.size() + " scores from " + q_score_range_file + ".");
			}

			q_score_range = new double[QR.size()][2];
			QRName = new String[QR.size()];

			for (int i = 0; i < QR.size(); i++) {
				ArrayList<String> qr = QR.get(i);
				QRName[i] = qr.get(0);
				q_score_range[i][0] = Double.parseDouble(qr.get(1));
				q_score_range[i][1] = Double.parseDouble(qr.get(2)); 
			}

			isQ = true;
		}

		if (par.MaCH_Dosage != null) {
			dosageFile = new String[1];
			dosageFile[0] = par.MaCH_Dosage;
			File f = new File(dosageFile[0]);
			if (!f.exists()) {
				System.err.println("could not open " + dosageFile[0] + ".");
				System.exit(0);
			}
			inforFile = new String[1];
			inforFile[0] = par.MaCH_Infor;
			f = new File(inforFile[0]);
			if (!f.exists()) {
				System.err.println("could not open " + dosageFile[0] + ".");
				System.exit(0);
			}

		} else {
			BufferedReader reader1 = FileProcessor.FileOpen(par.MaCH_Dosage_Batch);
			ArrayList<String> l1 = NewIt.newArrayList();
			String line = null;
			try {
				while((line=reader1.readLine())!=null) {
					if (line.length() == 0) continue;
					l1.add(line);
				}
			} catch (IOException e) {
				e.printStackTrace(System.err);
				System.exit(0);
			}
			dosageFile = (String[]) l1.toArray(new String[0]);
			for (int i = 0; i < dosageFile.length; i++) {
				File f = new File(dosageFile[i]);
				if (!f.exists()) {
					System.err.println("could not open "+ dosageFile[i] + ".");
					System.exit(0);
				}
			}

			BufferedReader reader2 = FileProcessor.FileOpen(par.MaCH_Infor_Batch);
			ArrayList<String> l2 = NewIt.newArrayList();
			try {
				while ((line=reader2.readLine())!=null) {
					if (line.length() == 0) continue;
					l2.add(line);
				}
			} catch (IOException e) {
				e.printStackTrace(System.err);
				System.exit(0);
			}
			inforFile = (String[]) l2.toArray(new String[0]);
			for (int i = 0; i < inforFile.length; i++) {
				File f = new File(inforFile[i]);
				if (!f.exists()) {
					System.err.println("could not open "+ inforFile[i] + ".");
					System.exit(0);
				}
			}
		}
	}

	public void makeProfile() {
		if (isQ) {
			multipleProfile();
		} else {
			singleProfile();
		}
	}

	public void multipleProfile() {

		int[] CC = new int[QRName.length];
		int[] CCSNP = new int[QRName.length];
		double[][] riskProfile= null;
		int ATGCLocus = 0;
		
		int sumSNPMapped = 0;
		
		for(int i = 0; i < inforFile.length; i++) {

			ArrayList<DosageInfor> SD = readDosageInfor(inforFile[i]);

			ArrayList<ArrayList<Double>> dosage = readDosage(dosageFile[i], i);

			double[][] rs = new double[dosage.size()][QRName.length];

			if (i == 0) riskProfile = new double[dosage.size()][QRName.length];

			int[] c = new int[QRName.length];
			int cSNP = 0;

			for (int j = 0; j < SD.size(); j++) {

				boolean qsL1Flag = false;
				boolean[] qsL2Flag = new boolean[QRName.length];
				Arrays.fill(qsL2Flag, false);
				
				DosageInfor di = SD.get(j);
				String snp = di.getSNP();
				String refA = di.getRefAllele();
				String refB = di.getSecAllele();

				if(di.isATGCLocus() && !par.keepATGCFlag) {
					ATGCLocus++;
					continue;
				}

				ScoreUnit su = null;
				if (Score.containsKey(snp) ) {
					su = Score.get(snp);
					if (su.isMissing()) {
						continue;
					}
					cSNP++;

					if(QS.containsKey(snp)) {
						QScore qs = QS.get(snp);
						if (!qs.isMissing()) {
							qsL1Flag = true;
							for (int k = 0; k < QRName.length; k++) {
								if (qs.getQScore() >= q_score_range[k][0] && qs.getQScore() <= q_score_range[k][1]) {
									qsL2Flag[k] = true;
									CCSNP[k]++;
								}
							}
						}
					}

				} else {
					continue;
				}
				
				if (!qsL1Flag) {
					continue;
				}


				for (int k = 0; k < dosage.size(); k++) {

					double locusScore = 0;
					if (refA.compareTo(su.getRefAllele()) == 0) {
						locusScore += dosage.get(k).get(j).doubleValue() * su.getScore(); 
					} else if (refB.compareTo(su.getRefAllele()) == 0) {
						locusScore += (2 - dosage.get(k).get(j).doubleValue()) * su.getScore();
					} else if (SNPMatch.Flip(refA).compareTo(su.getRefAllele()) == 0) {
						locusScore += dosage.get(k).get(j).doubleValue() * su.getScore();
					} else if (SNPMatch.Flip(refB).compareTo(su.getRefAllele()) == 0) {
						locusScore += (2 - dosage.get(k).get(j).doubleValue()) * su.getScore();
					} else {
						continue;
					}
						
					for (int l = 0; l < qsL2Flag.length; l++) {
						if (!qsL2Flag[l]) continue;
						rs[k][l] += locusScore;
						if (k == 0) {
							c[l]++;
							CC[l]++;
						}
					}
				}
			}

			sumSNPMapped += cSNP;
			System.out.println(dosageFile[i] + " mapped " + cSNP + " SNPs to the score file.");
			for (int j = 0; j < c.length; j++) {
				System.out.println("\t"+c[j] + " SNPs mapped to the range " + q_score_range[j][0] + " " + q_score_range[j][1]);
			}

			for (int j = 0; j < riskProfile.length; j++) {
				for (int k = 0; k < riskProfile[j].length; k++) {
					riskProfile[j][k] += rs[j][k];
				}
			}
		}

		for (int i = 0; i < riskProfile.length; i++) {
			for (int j = 0; j < riskProfile[i].length; j++) {
				if (CC[j] == 0) {
					riskProfile[i][j] = 0;
				} else {
					riskProfile[i][j] /= 2*CC[j];
				}
			}
		}

		if (!par.keepATGCFlag) {
			if (ATGCLocus > 1) {
				System.out.println(ATGCLocus + " ATGC loci were removed.");
			} else {
				System.out.println(ATGCLocus + " ATGC locus was removed.");
			}
		}
		
		System.out.println("In total " + sumSNPMapped + " SNPs mapped to the score file.");

		for (int i = 0; i < CCSNP.length; i++) {
			System.out.println(CCSNP[i] + " SNPs mapped to the range " + q_score_range[i][0] + " " + q_score_range[i][1]);
		}

		StringBuffer sbim = new StringBuffer();
		sbim.append(par.out);
		sbim.append(".profile");
		PrintStream predictorFile = FileProcessor.CreatePrintStream(sbim.toString());
		predictorFile.print("FID\tIID");
		for (int i = 0; i < QRName.length; i++) {
			predictorFile.print("\tSCORE"+ "." + QRName[i]);
		}
		predictorFile.println();

		for (int i = 0; i < riskProfile.length; i++) {
			String[] id = ID.get(i).split("->");
			predictorFile.print(id[0] + "\t" + id[1]);
			for (int j = 0; j < riskProfile[i].length; j++) {
				predictorFile.print("\t"+riskProfile[i][j]);
			}
			predictorFile.println();
		}
		predictorFile.close();
	}

	public void singleProfile() {

		int CC = 0;
		int CCSNP = 0;
		double[] riskProfile = null;

		int ATGCLocus = 0;
		for(int i = 0; i < inforFile.length; i++) {

			ArrayList<DosageInfor> SD = readDosageInfor(inforFile[i]);

			ArrayList<ArrayList<Double>> dosage = readDosage(dosageFile[i], i);

			double[] rs = new double[dosage.size()];
			
			if (i == 0) riskProfile = new double[dosage.size()];

			int c = 0;
			int cSNP = 0;
			for (int j = 0; j < SD.size(); j++) {
				DosageInfor di = SD.get(j);
				String snp = di.getSNP();
				String refA = di.getRefAllele();
				String refB = di.getSecAllele();
				
				if(di.isATGCLocus() && !par.keepATGCFlag) {
					ATGCLocus++;
					continue;
				}

				ScoreUnit su = null;
				if (Score.containsKey(snp) ) {
					su = Score.get(snp);
					if (su.isMissing()) {
						continue;
					}
					cSNP++;
					CCSNP++;
				} else {
					continue;
				}

				for (int k = 0; k < dosage.size(); k++) {
					if (refA.compareTo(su.getRefAllele()) == 0) {
						rs[k] += dosage.get(k).get(j).doubleValue() * su.getScore(); 
					} else if (refB.compareTo(su.getRefAllele()) == 0) {
						rs[k] += (2 - dosage.get(k).get(j).doubleValue()) * su.getScore();
					} else if (SNPMatch.Flip(refA).compareTo(su.getRefAllele()) == 0) {
						rs[k] += dosage.get(k).get(j).doubleValue() * su.getScore();
					} else if (SNPMatch.Flip(refB).compareTo(su.getRefAllele()) == 0) {
						rs[k] += (2 - dosage.get(k).get(j).doubleValue()) * su.getScore();
					} else {
						continue;
					}
					if (k == 0) {
						c++;
						CC++;
					}
				}
			}
			System.out.println(dosageFile[i] + " mapped " + cSNP + " SNPs to the score file, and " + c + " SNPs had scores.");

			for (int j = 0; j < rs.length; j++) {
				riskProfile[j] += rs[j];
			}
		}

		for (int i = 0; i < riskProfile.length; i++) {
			if (CC == 0) {
				riskProfile[i] = 0;
			} else {
				riskProfile[i] /= 2*CC;
			}
		}

		if (!par.keepATGCFlag) {
			if (ATGCLocus > 1) {
				System.out.println(ATGCLocus + " loci were removed.");
			} else {
				System.out.println(ATGCLocus + " locus was removed.");
			}
		}
		System.out.println("\nIn total " + CCSNP + " SNPs mapped to the score file, and " + CC + " SNPs had scores.");

		StringBuffer sbim = new StringBuffer();
		sbim.append(par.out);
		sbim.append(".profile");
		PrintStream predictorFile = FileProcessor.CreatePrintStream(sbim.toString());
		predictorFile.println("FID\tIID\tSCORE");
		for (int i = 0; i < riskProfile.length; i++) {
			String[] id = ID.get(i).split("->");
			predictorFile.println(id[0] + "\t" + id[1] + "\t"+riskProfile[i]);
		}
		predictorFile.close();
	}

	private ArrayList<DosageInfor> readDosageInfor(String file) {
		System.out.println("reading " + file);
		ArrayList<DosageInfor> SD = NewIt.newArrayList();

		BufferedReader readerDIFile = FileProcessor.FileOpen(file);
		String lineScore = null;

		try {
			lineScore = readerDIFile.readLine(); //header
			while((lineScore = readerDIFile.readLine())!=null) {
				DosageInfor di = new DosageInfor(lineScore);
				SD.add(di);
			}
		} catch (IOException e) {
			e.printStackTrace(System.err);
			System.exit(0);
		}
		return SD;
	}

	private ArrayList<ArrayList<Double>> readDosage(String file, int idx) {
		System.out.println("reading " + file);
		ArrayList<ArrayList<Double>> dosage = NewIt.newArrayList();
		BufferedReader readerDosageFile = FileProcessor.ZipFileOpen(file);
		String lineDosage = null;
		try {
			while((lineDosage = readerDosageFile.readLine())!=null) {
				String s[] = lineDosage.split(delim);
				if ( idx == 0) {
					ID.add(s[0]);
				}

				ArrayList<Double> d = NewIt.newArrayList();
				for (int i = 2; i < s.length; i++) {
					d.add(Double.parseDouble(s[i]));
				}
				dosage.add(d);
			}
		} catch (IOException e) {
			e.printStackTrace(System.err);
			System.exit(0);
		}

		return dosage;
	}
}
