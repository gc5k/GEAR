package profile;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import parameter.Parameter;
import profile.struct.QScore;
import profile.struct.ScoreUnit;
import test.Test;
import family.pedigree.file.SNP;
import family.plink.PLINKBinaryParser;
import family.plink.PLINKParser;
import family.popstat.GenotypeMatrix;
import family.qc.rowqc.SampleFilter;
import gear.util.FileProcessor;
import gear.util.NewIt;
import gear.util.SNPMatch;

public class RiskScore {
	private String delim = "\\s+";
	private GenotypeMatrix G1;

//	private int[][] comSNPIdx;
//	private double[][] allelefreq1;
//	private double[] N1;
//	private double[] N2;
//	private ArrayList<Boolean> flag;
	private HashMap<String, ScoreUnit> Score = NewIt.newHashMap();

	private String q_score_file;
	private String q_score_range_file;
	
	private String scoreFile;
//	private String[] title;
	private ArrayList<SNP> snpList1;
//	private ArrayList<Predictor> predictorList = NewIt.newArrayList();
	private HashMap<String, QScore> QS = NewIt.newHashMap();
	private double[][] q_score_range;
	private String[] QRName;
	private boolean isQ = false;

	ArrayList<Integer> scoreCoding = NewIt.newArrayList();

	private SampleFilter sf1;

	public RiskScore() {
		System.err.print("generating risk profile for genotypes.\n");
		
		initial();

		PLINKParser pp1 = null;
		if (Parameter.INSTANCE.hasBFileOption()) {
			pp1 = new PLINKBinaryParser (Parameter.INSTANCE.getBedFile(),
					                     Parameter.INSTANCE.getBimFile(),
					                     Parameter.INSTANCE.getFamFile());
		} else {
			System.err.println("did not specify files.");
			Test.LOG.append("did not specify files.\n");
			Test.printLog();
			System.exit(0);
		}
		pp1.Parse();

		sf1 = new SampleFilter(pp1.getPedigreeData(), pp1.getMapData());
		G1 = new GenotypeMatrix(sf1.getSample());
		snpList1 = sf1.getMapFile().getMarkerList();

	}

	private void initial() {

		// read score file
		scoreFile = Parameter.INSTANCE.scoreFile;
		BufferedReader readerScoreFile = FileProcessor.FileOpen(scoreFile);
		String lineScore = null;
		try {
			while ((lineScore = readerScoreFile.readLine()) != null) {
				if (lineScore.length() == 0)
					continue;
				ScoreUnit su = new ScoreUnit(lineScore);
				Score.put(su.getSNP(), su);
			}
		} catch (IOException e) {
			e.printStackTrace(System.err);
			System.exit(0);
		}

		System.err.println(Score.size() + " predictors are read from " + scoreFile + ".");

		// read q score file and q range file
		if (Parameter.INSTANCE.q_score_file != null && Parameter.INSTANCE.q_score_range_file != null) {

			System.out.println(Parameter.INSTANCE.q_score_file + " " + Parameter.INSTANCE.q_score_range_file);
			// q score file
			q_score_file = Parameter.INSTANCE.q_score_file;
			BufferedReader readerQScoreFile = FileProcessor
					.FileOpen(q_score_file);
			String lineQScore = null;
			try {
				while ((lineQScore = readerQScoreFile.readLine()) != null) {
					if (lineQScore.length() == 0)
						continue;
					lineQScore.trim();
					QScore qs = new QScore(lineQScore);
					QS.put(qs.getSNP(), qs);
				}
			} catch (IOException e) {
				e.printStackTrace(System.err);
				System.exit(0);
			}
			if (QS.size() == 0) {
				System.out.println("nothing has been selected in "
						+ q_score_file);
				System.exit(0);
			} else {
				System.out.println("read in " + QS.size() + " SNP scores from "
						+ q_score_file + ".");
			}

			// q range file
			q_score_range_file = Parameter.INSTANCE.q_score_range_file;
			BufferedReader readerQRangeFile = FileProcessor
					.FileOpen(q_score_range_file);
			String lineQRange = null;
			ArrayList<ArrayList<String>> QR = NewIt.newArrayList();
			try {
				while ((lineQRange = readerQRangeFile.readLine()) != null) {
					if (lineQRange.length() == 0)
						continue;
					String[] s = lineQRange.split(delim);
					if (s.length < 3)
						continue;
					ArrayList<String> qr = NewIt.newArrayList();
					qr.add(s[0]);
					qr.add(s[1]);
					qr.add(s[2]);
					QR.add(qr);
				}
			} catch (IOException e) {
				e.printStackTrace(System.err);
				System.exit(0);
			}
			if (QR.size() == 0) {
				System.out.println("nothing has been selected in "
						+ q_score_range_file);
			} else {
				System.out.println("read in " + QR.size() + " scores from "
						+ q_score_range_file + ".");
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

	}
	
	public void makeProfile() {
		if (isQ) {
			multipleProfile();
		} else {
			singleProfile();
		}
	}

	public void multipleProfile() {
		int Total = 0;
		int ATGCLocus = 0;
		int[] CCSNP = new int[QRName.length];
		int[][] GCInd = new int[G1.getGRow()][QRName.length];
		int[] matchScheme = new int[5];

		double[][] riskProfile = new double[G1.getGRow()][QRName.length];
//		int[] CC = new int[G1.getGRow()];
		for(int i = 0; i < snpList1.size(); i++) {

			boolean[] qsL2Flag = new boolean[QRName.length];
			Arrays.fill(qsL2Flag, false);

			SNP snp = snpList1.get(i);
			char a1 = snp.getRefAllele();
			char a2 = snp.getSecAllele();
			boolean isATGC = SNPMatch.Confusion(a1, a2);

			ScoreUnit su = null;
			boolean isMatch = true;
			double sc = 0;
			if (Score.containsKey(snp.getName()) ) {
				Total++;
				if (isATGC) {
					ATGCLocus++;
					if (!Parameter.INSTANCE.keepATGC()) {
						continue;
					}
				}
				su = Score.get(snp.getName());
				if (su.isMissing()) {
					continue;
				}

				if(QS.containsKey(snp.getName())) {
					QScore qs = QS.get(snp.getName());
					if (!qs.isMissing()) {
						for (int k = 0; k < QRName.length; k++) {
							if (qs.getQScore() >= q_score_range[k][0] && qs.getQScore() <= q_score_range[k][1]) {
								qsL2Flag[k] = true;
								CCSNP[k]++;
							}
						}
					}
				} else {
					continue;
				}

				if (su.getRefAllele().compareTo(Character.toString(a1)) == 0) {
					isMatch = true;
					matchScheme[0]++;

				} else if (su.getRefAllele().compareTo(Character.toString(a2)) == 0) {
					isMatch = false;
					matchScheme[1]++;

				} else if (su.getRefAllele().compareTo(SNPMatch.Flip(Character.toString(a1)))==0){
					isMatch = true;
					matchScheme[2]++;

				} else if (su.getRefAllele().compareTo(SNPMatch.Flip(Character.toString(a2)))==0) {
					isMatch = false;
					matchScheme[3]++;

				} else {
					matchScheme[4]++;
					continue;
				}

				sc = su.getScore();
				if(Parameter.INSTANCE.getTranFunction() == parameter.RegressionModel.LOGIT) {
					if(isMatch) {
						sc = Math.log(sc);
					} else {
						sc = -1 * Math.log(sc);
					}
				} else {
					if(!isMatch) {
						sc = -1 * sc; 
					}
				}
			} else {// this snp is not in the predictor panel;
				continue;
			}

			for(int k = 0; k < qsL2Flag.length; k++) {
				if (!qsL2Flag[k]) continue;
				CCSNP[k]++;
				for(int j = 0; j < G1.getGRow(); j++) {
					if (G1.getAdditiveScore(j, i)!=GenotypeMatrix.missing) {
						riskProfile[j][k] += sc * (2-G1.getAdditiveScore(j, i));
						GCInd[j][k]++;
					}
				}
			}
		}

		for (int i = 0; i < riskProfile.length; i++) {
			for (int j = 0; j < QRName.length; j++) {
				if (GCInd[i][j] == 0) {
					riskProfile[i][j] = 0;
				} else {
					riskProfile[i][j] /= 2*GCInd[i][j];
				}
			}
		}

		System.out.println(Total + " SNPs were mapped to the score file.");
		if (ATGCLocus > 1) {
			if(Parameter.INSTANCE.keepATGC()) {
				System.out.println(ATGCLocus + " ATGC loci were detected.");
			} else {
				System.out.println(ATGCLocus + " ATGC loci were removed.");
			}
		} else {
			if(Parameter.INSTANCE.keepATGC()) {
				System.out.println(ATGCLocus + " ATGC locus was detected.");
			} else {
				System.out.println(ATGCLocus + " ATGC locus was removed.");
			}
		}
		for (int i = 0; i < QRName.length; i++) {
			System.out.println(CCSNP[i] + " SNPs have mapped into " + QRName[i] + ".");
		}

		for (int i = 0; i < 4; i++) {
			System.out.println(matchScheme[i] + " SNPs match Scheme " + (1+i));
		}

		StringBuffer sbim = new StringBuffer();
		sbim.append(Parameter.INSTANCE.out);
		sbim.append(".profile");
		PrintStream predictorFile = FileProcessor.CreatePrintStream(sbim.toString());
		
		predictorFile.print("FID\tIID\tPHENO");
		for (int i = 0; i < QRName.length; i++) {
			predictorFile.print("\tScore."+ QRName[i]);
		}
		predictorFile.println();
		for (int i = 0; i < riskProfile.length; i++) {
			predictorFile.print(sf1.getSample().get(i).getFamilyID() + "\t" + sf1.getSample().get(i).getIndividualID() + "\t"+ sf1.getHukouBook().get(i).getCol6());
			for (int j = 0; j < riskProfile[i].length; j++) {
				predictorFile.print("\t" + riskProfile[i][j]);
			}
			predictorFile.println();
		}
		predictorFile.close();
	}

	public void singleProfile() {

		int Total = 0;
		int CCSNP = 0;
		int ATGCLocus = 0;
		int[] matchScheme = new int[5];

		ArrayList<ArrayList<String>> s4 = NewIt.newArrayList();
		for(int i = 0; i < 4; i++) {
			ArrayList<String> s = NewIt.newArrayList();
			s4.add(s);
		}
		double[] riskProfile = new double[G1.getGRow()];
		int[] GCInd = new int[G1.getGRow()];
//		int[] CC = new int[G1.getGRow()];
		for(int i = 0; i < snpList1.size(); i++) {
			SNP snp = snpList1.get(i);
			char a1 = snp.getRefAllele();
			char a2 = snp.getSecAllele();
			boolean isATGC = SNPMatch.Confusion(a1, a2);

			ScoreUnit su = null;
			boolean isMatch = true;
			double sc = 0;
			if (Score.containsKey(snp.getName()) ) {
				Total++;
				if (isATGC) {
					ATGCLocus++;
					if (!Parameter.INSTANCE.keepATGC()) {
						continue;
					}
				}
				su = Score.get(snp.getName());
//				System.err.println(snp.getName()+ "\t" + Score.get(snp.getName()).getSNP());
				if (su.isMissing()) {
					continue;
				}

				if (su.getRefAllele().compareTo(Character.toString(a1)) == 0) {
					isMatch = true;
					matchScheme[0]++;
					ArrayList<String> s = s4.get(0);
					s.add(snp.getName());
				} else if (su.getRefAllele().compareTo(Character.toString(a2)) == 0) {
					isMatch = false;
					matchScheme[1]++;
					ArrayList<String> s = s4.get(1);
					s.add(snp.getName());

				} else if (su.getRefAllele().compareTo(SNPMatch.Flip(Character.toString(a1)))==0){
					isMatch = true;
					matchScheme[2]++;
					ArrayList<String> s = s4.get(2);
					s.add(snp.getName());

				} else if (su.getRefAllele().compareTo(SNPMatch.Flip(Character.toString(a2)))==0) {
					isMatch = false;
					matchScheme[3]++;
					ArrayList<String> s = s4.get(3);
					s.add(snp.getName());

				} else {
					matchScheme[4]++;
					continue;
				}
				CCSNP++;

				sc = su.getScore();
				if (Parameter.INSTANCE.getTranFunction() == parameter.RegressionModel.LOGIT) {//logit s
					if(isMatch) {
						sc = Math.log(sc);
					} else {
						sc = -1 * Math.log(sc);
					}
				}
			} else {// this snp is not in the predictor panel;
				continue;
			}

			for(int j = 0; j < G1.getGRow(); j++) {
				if(G1.getAdditiveScore(j, i) != GenotypeMatrix.missing) {
					riskProfile[j] += sc * (2-G1.getAdditiveScore(j, i));
					GCInd[j]++;
				}
			}
		}

		for (int i = 0; i < riskProfile.length; i++) {
			if (GCInd[i] == 0) {
				riskProfile[i] = 0;
			} else {
				riskProfile[i] /= 2*GCInd[i];
			}
		}

		System.out.println(Total + " SNPs were mapped to the score file.");
		if (ATGCLocus > 1) {
			if (Parameter.INSTANCE.keepATGC()) {
				System.out.println(ATGCLocus + " ATGC loci were detected.");
			} else {
				System.out.println(ATGCLocus + " ATGC loci were removed.");
			}
		} else {
			if (Parameter.INSTANCE.keepATGC()) {
				System.out.println(ATGCLocus + " ATGC Locus were detected.");
			} else {
				System.out.println(ATGCLocus + " ATGC locus was removed.");				
			}
		}
		System.out.println(CCSNP + " SNPs have score in the score file.");

		for (int i = 0; i < 4; i++) {
			System.out.println(matchScheme[i] + " SNPs match Scheme " + (1+i));
		}

		StringBuffer sbim = new StringBuffer();
		sbim.append(Parameter.INSTANCE.out);
		sbim.append(".profile");
		PrintStream predictorFile = FileProcessor.CreatePrintStream(sbim.toString());
		predictorFile.println("FID\tIID\tPHENO\tSCORE");
		for (int i = 0; i < riskProfile.length; i++) {
			predictorFile.println(sf1.getSample().get(i).getFamilyID() + "\t" + sf1.getSample().get(i).getIndividualID() + "\t"+ sf1.getHukouBook().get(i).getCol6() + "\t" + riskProfile[i]);
		}
		predictorFile.close();
/*		
		PrintStream pf = FileProcessor.CreatePrintStream("s1.txt");
		ArrayList<String> S = s4.get(0);
		for (int i = 0; i < S.size(); i++) {
			pf.println(S.get(i));
		}
		pf.close();
		
		pf = FileProcessor.CreatePrintStream("s2.txt");
		S = s4.get(1);
		for (int i = 0; i < S.size(); i++) {
			pf.println(S.get(i));
		}
		pf.close();

		pf = FileProcessor.CreatePrintStream("s3.txt");
		S = s4.get(2);
		for (int i = 0; i < S.size(); i++) {
			pf.println(S.get(i));
		}
		pf.close();
		
		pf = FileProcessor.CreatePrintStream("s4.txt");
		S = s4.get(3);
		for (int i = 0; i < S.size(); i++) {
			pf.println(S.get(i));
		}
		pf.close();
*/
	}
}
