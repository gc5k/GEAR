package gear.subcommands.wgrmA;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import gear.ConstValues;
import gear.data.Person;
import gear.data.SubjectID;
import gear.family.GenoMatrix.GenotypeMatrix;
import gear.family.pedigree.PersonIndex;
import gear.family.plink.PLINKParser;
import gear.qc.sampleqc.SampleFilter;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.BufferedReader;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;
import gear.util.pop.PopStat;
import gear.util.stat.CommonMath;

public class WGRMACommandImpl extends CommandImpl {
	private HashMap<String, Float> scores = NewIt.newHashMap(); // name-to-score

	private GenotypeMatrix pGM;

	private float[][] allelefreq;
	private float[] allelevar;
	private float[] weight;

	private WGRMACommandArguments wgrmArgs;
	private SampleFilter sf;
	private HashSet<SubjectID> subjectSet =	NewIt.newHashSet();

	private int[][] Gcnt = null;
	private int grmTriangleSize;
	private int row;
	private int col;

	ArrayList<ArrayList<Integer>> missList = NewIt.newArrayList();

	@Override
	public void execute(CommandArguments cmdArgs) {
		wgrmArgs = (WGRMACommandArguments) cmdArgs;

		PLINKParser pp = PLINKParser.parse(wgrmArgs);
		sf = new SampleFilter(pp.getPedigreeData(), cmdArgs);
		pGM = new GenotypeMatrix(sf.getSample(), pp.getMapData(), cmdArgs);

		Logger.printUserLog("");
		Logger.printUserLog("Calculating allele frequencies...");
		allelefreq = PopStat.calAlleleFrequencyFloat(pGM);
		allelevar = PopStat.calGenoVarianceFloat(pGM);
		weight = new float[pGM.getNumMarker()];

		if(wgrmArgs.isWeight() || wgrmArgs.isVanRaden()) {
			prepareWeight();
		} else {
			Arrays.fill(weight, 1);
		}

		if (wgrmArgs.isInbredList()) {
			readInbredList();
		}

		Gcnt = new int[pGM.getNumIndivdial()][pGM.getNumIndivdial()];

		makeGA();

		if (wgrmArgs.isDom()) {
			makeGD();			
		}
//		makeAddScore();
//
//		if (grmArgs.isDom()) {
//			makeDomScore();
//		}
	}

	private void prepareWeight() {
		Logger.printUserLog("Preparing weights...");

		if (wgrmArgs.isWeight()) {
			BufferedReader reader = BufferedReader.openTextFile(wgrmArgs.getWeightFile(), "weight");

			String[] tokens = reader.readTokens(2);
			int cnt = 0;
			while (tokens != null) {
				Double.parseDouble(tokens[1]);
				scores.put(tokens[0], Float.parseFloat(tokens[1]));
				cnt++;
				tokens = reader.readTokens(2);
			}
			reader.close();
			Logger.printUserLog(cnt + " weights have been read.");

			int scnt = 0;
			for (int i = 0; i < pGM.getSNPList().size(); i++) {
				float s = 0;
				if (scores.containsKey(pGM.getSNPList().get(i).getName())) {
					s = scores.get(pGM.getSNPList().get(i).getName());
					scnt++;
				}
				weight[i] = s;
			}
			Logger.printUserLog(scnt + " weights have been mapped.");
		} else if (wgrmArgs.isVanRaden()) {
			for (int i = 0; i < pGM.getSNPList().size(); i++) {

				if (allelefreq[i][1] == Double.NaN || allelefreq[i][0] == 0 || allelefreq[i][1] == 0 || allelevar[i] == 0) continue;

				weight[i] = (float) (wgrmArgs.isInbred()? Math.sqrt(4 * allelefreq[i][1] * (1 - allelefreq[i][1])):Math.sqrt(2 * allelefreq[i][1] * (1 - allelefreq[i][1])));
				if (wgrmArgs.isAdjVar()) weight[i] = (float) Math.sqrt(allelevar[i]);
			}
		}
	}

	private void makeGA() {
		Logger.printUserLog("Making additive genetic relatedness matrix...");
		long startNanoTime = System.nanoTime();
		float[][] gMat = new float[pGM.getNumIndivdial()][pGM.getNumMarker()];
		for (int i = 0; i < pGM.getNumIndivdial(); i++) {
			ArrayList<Integer> mL = NewIt.newArrayList();
			missList.add(mL);
		}

		Logger.printUserLog("Standardizing genotypes...");
		float W = 0;
		for (int i = 0; i < pGM.getNumMarker(); i++) {
			if (allelefreq[i][1] == Float.NaN || allelefreq[i][0] == 0 || allelefreq[i][1] == 0 || allelevar[i] == 0) continue;
			W += weight[i]*weight[i];

			float f = allelefreq[i][1];
			float de = 0;
			if (wgrmArgs.isInbred()) {
				de = (float) Math.sqrt(4*allelefreq[i][1]*allelefreq[i][0]);
			} else if (wgrmArgs.isAdjVar()) {
				de = (float) Math.sqrt(allelevar[i]);
			} else {
				de = (float) Math.sqrt(2*allelefreq[i][1]*allelefreq[i][0]);
			}

			for (int j = 0; j < pGM.getNumIndivdial(); j++) {
				int g = pGM.getAdditiveScore(j,i);
				if (g == Person.MissingGenotypeCode) {
					ArrayList<Integer> mL = missList.get(j);
					mL.add(i);
					missList.set(j, mL);
					gMat[j][i] = 0;
				} else {
					gMat[j][i] = (float) ((g*1.0-2*f)/de * weight[i]);
				}
			}
		}

		grmTriangleSize = CommonMath.calculateTriangleSize(gMat.length);
		if (grmTriangleSize < 15) {
			Logger.printUserError("Too small sample size. GEAR quit.");
			System.exit(1);
		}

		float[][] GA = new float[gMat.length][gMat.length];
		
		ProgressGA progress = new ProgressGA();
		progress.start();

		for (row = 0; row < gMat.length; row++) {
			float ws1 = W;
			ArrayList<Integer> mL1 = missList.get(row);
			for(int j = 0; j < mL1.size(); j++) {
				ws1 -= weight[mL1.get(j)]*weight[mL1.get(j)];		
			}

			for (col = 0; col <= row; col++) {
				float ws2 = ws1;
				ArrayList<Integer> mL2 = missList.get(col);

				for(int k = 0; k < mL2.size(); k++) {
					ws2 -= weight[mL2.get(k)]*weight[mL2.get(k)];
				}

				int mLoci = gMat[row].length;
				int cnt = 0;
				if ( mL1.size() <= mL2.size() ) {
					for (int k = 0; k < mL1.size(); k++) {
						if (Collections.binarySearch(mL2, mL1.get(k)) >=0) {
							ws2 += weight[mL1.get(k)]*weight[mL1.get(k)];
							cnt++;
						}
					}
				} else {
					for (int k = 0; k < mL2.size(); k++) {
						if (Collections.binarySearch(mL1, mL2.get(k)) >=0) {
							ws2 += weight[mL2.get(k)]*weight[mL2.get(k)];
							cnt++;
						}
					}
				}

				mLoci -= (mL1.size() + mL2.size() - cnt);
				
				float s = 0; 
				for (int l = 0; l < gMat[row].length; l++) {
					s += gMat[row][l] * gMat[col][l];
				}
				GA[row][col] = s/ws2;
				Gcnt[row][col] = mLoci;
			}
		}

		try {
			progress.join();
		} catch (InterruptedException e) { }
		Logger.printUserLog("");

		double grmMean = 0;
		double grmSq = 0;

		StringBuffer sb = new StringBuffer();
		sb.append(wgrmArgs.getOutRoot());
		PrintStream grm = null;
		BufferedWriter grmGZ = null;
		if (wgrmArgs.isGZ()) {
			sb.append(".grm.gz");
			grmGZ = FileUtil.ZipFileWriter(sb.toString());
		} else {
			sb.append(".grm.txt");
			grm = FileUtil.CreatePrintStream(sb.toString());
		}

		int cnt = 0;
		for (int i = 0; i < Gcnt.length; i++) {
			for (int j = 0; j <= i; j++) {
				if (i != j) {
					grmMean += GA[i][j];
					grmSq += GA[i][j] * GA[i][j];
					cnt++;
				}
				if (wgrmArgs.isGZ()) {
					try {
						grmGZ.append((i + 1) + "\t" + (j + 1) + "\t" + Gcnt[i][j] + "\t" + GA[i][j] + "\n");
					} catch (IOException e) {
						Logger.handleException(e,
								"error in writing '" + sb.toString() + "' for " + (i + 1) + " " + (j + 1) + ".");
					}
				} else {
					grm.println((i + 1) + "\t" + (j + 1) + "\t" + Gcnt[i][j] + "\t" + GA[i][j]);
				}
			}
		}

		if (wgrmArgs.isGZ()) {
			try {
				grmGZ.close();
			} catch (IOException e) {
				Logger.handleException(e, " error in closing '" + sb.toString() + "'.");
			}
		} else {
			grm.close();
		}
		Logger.printUserLog("Writing GRM scores into '" + sb.toString() + "'.");
		StringBuffer sb_id = new StringBuffer();
		sb_id.append(wgrmArgs.getOutRoot());
		PrintStream grm_id = null;
		sb_id.append(".grm.id");
		grm_id = FileUtil.CreatePrintStream(sb_id.toString());

		ArrayList<PersonIndex> PI = pGM.getSample();
		for (int i = 0; i < PI.size(); i++) {
			grm_id.println(PI.get(i).getFamilyID() + "\t" + PI.get(i).getIndividualID());
		}
		grm_id.close();
		Logger.printUserLog("Writing " + PI.size() + " individuals' information into '" + sb_id.toString() + "'.");

		grmMean /= cnt;
		grmSq /= cnt;
		double Effective_sample = -1 / grmMean + 1;
		double grmSD = (grmSq - grmMean * grmMean) * cnt / (cnt - 1);
		double Effective_marker = 1 / grmSD;

		DecimalFormat df = new DecimalFormat("0.0000");
		DecimalFormat dfE = new DecimalFormat("0.00E0");
		if (Math.abs(grmMean) > 0.0001) {
			Logger.printUserLog("Mean of genetic relatedness is: " + df.format(grmMean));
		} else {
			Logger.printUserLog("Mean of genetic relatedness is: " + dfE.format(grmMean));
		}

		if (Math.abs(grmSD) > 0.0001) {
			Logger.printUserLog("Sampling variance of genetic relatedness is: " + df.format(grmSD));
		} else {
			Logger.printUserLog("Sampling variance of genetic relatedness is: " + dfE.format(grmSD));
		}

		if (Math.abs(Effective_marker) > 0.0001) {
			Logger.printUserLog("Effective sample size is: " + df.format(Effective_sample));
		} else {
			Logger.printUserLog("Effective sample size is: " + dfE.format(Effective_sample));
		}

		if (Math.abs(Effective_marker) > 0.0001) {
			Logger.printUserLog("Effective number of genome segments is: " + df.format(Effective_marker));
		} else {
			Logger.printUserLog("Effective number of genome segments is: " + dfE.format(Effective_marker));
		}

		if(wgrmArgs.isGUI()) {
			PrintStream gui_file = null;
			gui_file = FileUtil.CreatePrintStream(wgrmArgs.getOutRoot()+".grm.gui");
			gui_file.println(df.format(Effective_sample) + "\t" + PI.size());
			gui_file.println(df.format(Effective_marker) + "\t" + allelefreq.length);
			gui_file.close();
		}
		
		Logger.printElapsedTime(startNanoTime, "make additive genetic relatedness matrix");
	}
	
	class ProgressGA extends Thread {
		public void run() {
			int calculatedElements;
			do {
				calculatedElements = CommonMath.calculateTriangleSize(row) + col + 1;
				float percentage = Math.min(100f, (float)calculatedElements / grmTriangleSize * 100f);
				System.out.print(String.format(
						"[INFO] Calculating additive genetic relatedness matrix, %.2f%% completed...\r", percentage));
				try {
					Thread.sleep(1000);
				} catch (InterruptedException  e) { }
			} while (calculatedElements < grmTriangleSize);
		}
	}

	private void makeGD() {
		float[][] gMat = new float[pGM.getNumIndivdial()][pGM.getNumMarker()];

		for (int i = 0; i < pGM.getNumMarker(); i++) {
			if (allelefreq[i][1] == Float.NaN || allelefreq[i][0] == 0 || allelefreq[i][1] == 0 || allelevar[i] == 0) continue;

			float f = allelefreq[i][1];

			for (int j = 0; j < pGM.getNumIndivdial(); j++) {
				int g = pGM.getAdditiveScore(j,i);
				if (g == Person.MissingGenotypeCode) {
					ArrayList<Integer> mL = missList.get(j);
					mL.add(i);
					missList.set(j, mL);
					gMat[j][i] = 0;
				} else {					
					if (g == 0) {
						gMat[j][i] = (0 - 2*f*f) / (2*f*(1-f));
					} else if (g == 1) {
						gMat[j][i] = (2*f - 2*f*f) /(2*f*(1-f));
					} else {
						gMat[j][i] = ((4*f - 2) -(2*f*f))/ (2*f*(1-f));
					}
				}
			}
		}

		float[][] GD = new float[gMat.length][gMat.length];

		for (int i = 0; i < gMat.length; i++) {
			for (int j = 0; j <= i; j++) {
				float s = 0;
				for (int l = 0; l < gMat[i].length; l++) {
					s += gMat[i][l] * gMat[j][l];
				}
				GD[i][j] = s/Gcnt[i][j];
			}
		}

		double grmDomMean = 0;
		double grmDomSq = 0;

		StringBuffer sb = new StringBuffer();
		sb.append(wgrmArgs.getOutRoot());
		PrintStream grm = null;
		BufferedWriter grmGZ = null;
		if (wgrmArgs.isGZ()) {
			sb.append(".dom.grm.gz");
			grmGZ = FileUtil.ZipFileWriter(sb.toString());
		} else {
			sb.append(".dom.grm.txt");
			grm = FileUtil.CreatePrintStream(sb.toString());
		}

		int cnt = 0;
		for (int i = 0; i < GD.length; i++) {
			for (int j = 0; j <= i; j++) {
				if (i != j) {
					grmDomMean += GD[i][j];
					grmDomSq += GD[i][j] * GD[i][j];
					cnt++;
				}
				if (wgrmArgs.isGZ()) {
					try {
						grmGZ.append((i + 1) + "\t" + (j + 1) + "\t" + Gcnt[i][j] + "\t" + GD[i][j] + "\n");
					} catch (IOException e) {
						Logger.handleException(e,
								"error in writing '" + sb.toString() + "' for " + (i + 1) + " " + (j + 1) + ".");
					}
				} else {
					grm.println((i + 1) + "\t" + (j + 1) + "\t" + Gcnt[i][j] + "\t" + GD[i][j]);
				}
			}
		}

		if (wgrmArgs.isGZ()) {
			try {
				grmGZ.close();
			} catch (IOException e) {
				Logger.handleException(e, " error in closing '" + sb.toString() + "'.");
			}
		} else {
			grm.close();
		}
		Logger.printUserLog("Writing GRM dominance scores into '" + sb.toString() + "'.");
		StringBuffer sb_id = new StringBuffer();
		sb_id.append(wgrmArgs.getOutRoot());
		PrintStream grm_id = null;
		sb_id.append(".dom.grm.id");
		grm_id = FileUtil.CreatePrintStream(sb_id.toString());

		ArrayList<PersonIndex> PI = sf.getSample();
		for (int i = 0; i < PI.size(); i++) {
			grm_id.println(PI.get(i).getFamilyID() + "\t" + PI.get(i).getIndividualID());
		}
		grm_id.close();
		Logger.printUserLog("Writing " + PI.size() + " individuals' information into '" + sb_id.toString() + "'.");

		grmDomMean /= cnt;
		grmDomSq /= cnt;
		double Effective_DomSample = -1 / grmDomMean;
		double grmDomSD = (grmDomSq - grmDomMean * grmDomMean) * cnt / (cnt - 1);
		double Effeictive_DomMarker = 1 / grmDomSD;

		DecimalFormat df = new DecimalFormat("0.0000");
		DecimalFormat dfE = new DecimalFormat("0.00E0");
		if (Math.abs(grmDomMean) > 0.0001) {
			Logger.printUserLog("Mean of dominance genetic relatedness is: " + df.format(grmDomMean));
		} else {
			Logger.printUserLog("Mean of dominance genetic relatedness is: " + dfE.format(grmDomMean));
		}

		if (Math.abs(grmDomSD) > 0.0001) {
			Logger.printUserLog("Sampling variance of dominance genetic relatedness is: " + df.format(grmDomSD));
		} else {
			Logger.printUserLog("Sampling variance of dominance genetic relatedness is: " + dfE.format(grmDomSD));
		}

		if (Math.abs(Effeictive_DomMarker) > 0.0001) {
			Logger.printUserLog("Effective dominance sample size is: " + df.format(Effective_DomSample));
		} else {
			Logger.printUserLog("Effective dominance sample size is: " + dfE.format(Effective_DomSample));
		}

		if (Math.abs(Effeictive_DomMarker) > 0.0001) {
			Logger.printUserLog("Effective number of dominance genome segments is: " + df.format(Effeictive_DomMarker));
		} else {
			Logger.printUserLog(
					"Effective number of dominance genome segments is: " + dfE.format(Effeictive_DomMarker));
		}
	}

	private void readInbredList() {
		BufferedReader reader = BufferedReader.openTextFile(wgrmArgs.getInbredList(), "inbred list");
		int numCols = 2;

		String[] tokens = reader.readTokensAtLeast(numCols);

		if (tokens == null)	{
			Logger.printUserError("The file '" + wgrmArgs.getInbredList() + "' is empty.");
			System.exit(1);
		}

		do {
			SubjectID subjectID = new SubjectID(/*famID*/tokens[0], /*indID*/tokens[1]);
			if (subjectSet.contains(subjectID)) {
				Logger.printUserLog("Subject " + subjectID + " is duplicated.");
			} else {
				subjectSet.add(subjectID);
			}
		} while ((tokens = reader.readTokensAtLeast(numCols)) != null);
		reader.close();

		ArrayList<PersonIndex> sList= pGM.getSample();
		int cCnt = 0;
		for(int i = 0; i < sList.size(); i++) {
			SubjectID sID = new SubjectID(sList.get(i).getFamilyID(), sList.get(i).getIndividualID());
			if (subjectSet.contains(sID)) cCnt++;
		}
		Logger.printUserLog("Read " + subjectSet.size() + " individuals from '" + wgrmArgs.getInbredList() + "'.");
		Logger.printUserLog("Matched '" + cCnt + "' individuals.");

		if (cCnt < ConstValues.TooFewSample) {
			Logger.printUserLog("Too few individuals left for analysis <" + ConstValues.TooFewSample+". GEAR quit.");
			System.exit(1);
		}
	}

	public void makeDomScore() {
		double grmDomMean = 0;
		double grmDomSq = 0;

		StringBuffer sb = new StringBuffer();
		sb.append(wgrmArgs.getOutRoot());
		PrintStream grm = null;
		BufferedWriter grmGZ = null;
		if (wgrmArgs.isGZ()) {
			sb.append(".dom.grm.gz");
			grmGZ = FileUtil.ZipFileWriter(sb.toString());
		} else {
			sb.append(".dom.grm.txt");
			grm = FileUtil.CreatePrintStream(sb.toString());
		}

		int cnt = 0;
		for (int i = 0; i < pGM.getGRow(); i++) {
			for (int j = 0; j <= i; j++) {
				double[] d = GRMDomScore(i, j);
				if (i != j) {
					grmDomMean += d[1];
					grmDomSq += d[1] * d[1];
					cnt++;
				}
				if (wgrmArgs.isGZ()) {
					try {
						grmGZ.append((i + 1) + "\t" + (j + 1) + "\t" + d[0] + "\t" + d[1] + "\n");
					} catch (IOException e) {
						Logger.handleException(e,
								"error in writing '" + sb.toString() + "' for " + (i + 1) + " " + (j + 1) + ".");
					}
				} else {
					grm.println((i + 1) + "\t" + (j + 1) + "\t" + d[0] + "\t" + d[1]);
				}
			}
		}

		if (wgrmArgs.isGZ()) {
			try {
				grmGZ.close();
			} catch (IOException e) {
				Logger.handleException(e, " error in closing '" + sb.toString() + "'.");
			}
		} else {
			grm.close();
		}
		Logger.printUserLog("Writing GRM dominance scores into '" + sb.toString() + "'.");
		StringBuffer sb_id = new StringBuffer();
		sb_id.append(wgrmArgs.getOutRoot());
		PrintStream grm_id = null;
		sb_id.append(".dom.grm.id");
		grm_id = FileUtil.CreatePrintStream(sb_id.toString());

		ArrayList<PersonIndex> PI = sf.getSample();
		for (int i = 0; i < PI.size(); i++) {
			grm_id.println(PI.get(i).getFamilyID() + "\t" + PI.get(i).getIndividualID());
		}
		grm_id.close();
		Logger.printUserLog("Writing " + PI.size() + " individuals' information into '" + sb_id.toString() + "'.");

		grmDomMean /= cnt;
		grmDomSq /= cnt;
		double Effective_DomSample = -1 / grmDomMean;
		double grmDomSD = (grmDomSq - grmDomMean * grmDomMean) * cnt / (cnt - 1);
		double Effeictive_DomMarker = 1 / grmDomSD;

		DecimalFormat df = new DecimalFormat("0.0000");
		DecimalFormat dfE = new DecimalFormat("0.00E0");
		if (Math.abs(grmDomMean) > 0.0001) {
			Logger.printUserLog("Mean of dominance genetic relatedness is: " + df.format(grmDomMean));
		} else {
			Logger.printUserLog("Mean of dominance genetic relatedness is: " + dfE.format(grmDomMean));
		}

		if (Math.abs(grmDomSD) > 0.0001) {
			Logger.printUserLog("Sampling variance of dominance genetic relatedness is: " + df.format(grmDomSD));
		} else {
			Logger.printUserLog("Sampling variance of dominance genetic relatedness is: " + dfE.format(grmDomSD));
		}

		if (Math.abs(Effeictive_DomMarker) > 0.0001) {
			Logger.printUserLog("Effective dominance sample size is: " + df.format(Effective_DomSample));
		} else {
			Logger.printUserLog("Effective dominance sample size is: " + dfE.format(Effective_DomSample));
		}

		if (Math.abs(Effeictive_DomMarker) > 0.0001) {
			Logger.printUserLog("Effective number of dominance genome segments is: " + df.format(Effeictive_DomMarker));
		} else {
			Logger.printUserLog(
					"Effective number of dominance genome segments is: " + dfE.format(Effeictive_DomMarker));
		}
	}


	public void makeAddScore() {
		double grmMean = 0;
		double grmSq = 0;

		StringBuffer sb = new StringBuffer();
		sb.append(wgrmArgs.getOutRoot());
		PrintStream grm = null;
		BufferedWriter grmGZ = null;
		if (wgrmArgs.isGZ()) {
			sb.append(".grm.gz");
			grmGZ = FileUtil.ZipFileWriter(sb.toString());
		} else {
			sb.append(".grm.txt");
			grm = FileUtil.CreatePrintStream(sb.toString());
		}

		int cnt = 0;
		for (int i = 0; i < pGM.getGRow(); i++) {
			for (int j = 0; j <= i; j++) {
				double[] s = GRMScore(i, j);
				if (i != j) {
					grmMean += s[1];
					grmSq += s[1] * s[1];
					cnt++;
				}
				if (wgrmArgs.isGZ()) {
					try {
						grmGZ.append((i + 1) + "\t" + (j + 1) + "\t" + s[0] + "\t" + s[1] + "\n");
					} catch (IOException e) {
						Logger.handleException(e,
								"error in writing '" + sb.toString() + "' for " + (i + 1) + " " + (j + 1) + ".");
					}
				} else {
					grm.println((i + 1) + "\t" + (j + 1) + "\t" + s[0] + "\t" + s[1]);
				}
			}
		}

		if (wgrmArgs.isGZ()) {
			try {
				grmGZ.close();
			} catch (IOException e) {
				Logger.handleException(e, " error in closing '" + sb.toString() + "'.");
			}
		} else {
			grm.close();
		}
		Logger.printUserLog("Writing GRM scores into '" + sb.toString() + "'.");
		StringBuffer sb_id = new StringBuffer();
		sb_id.append(wgrmArgs.getOutRoot());
		PrintStream grm_id = null;
		sb_id.append(".grm.id");
		grm_id = FileUtil.CreatePrintStream(sb_id.toString());

		ArrayList<PersonIndex> PI = pGM.getSample();
		for (int i = 0; i < PI.size(); i++) {
			grm_id.println(PI.get(i).getFamilyID() + "\t" + PI.get(i).getIndividualID());
		}
		grm_id.close();
		Logger.printUserLog("Writing " + PI.size() + " individuals' information into '" + sb_id.toString() + "'.");

		grmMean /= cnt;
		grmSq /= cnt;
		double Effective_sample = -1 / grmMean + 1;
		double grmSD = (grmSq - grmMean * grmMean) * cnt / (cnt - 1);
		double Effective_marker = 1 / grmSD;

		DecimalFormat df = new DecimalFormat("0.0000");
		DecimalFormat dfE = new DecimalFormat("0.00E0");
		if (Math.abs(grmMean) > 0.0001) {
			Logger.printUserLog("Mean of genetic relatedness is: " + df.format(grmMean));
		} else {
			Logger.printUserLog("Mean of genetic relatedness is: " + dfE.format(grmMean));
		}

		if (Math.abs(grmSD) > 0.0001) {
			Logger.printUserLog("Sampling variance of genetic relatedness is: " + df.format(grmSD));
		} else {
			Logger.printUserLog("Sampling variance of genetic relatedness is: " + dfE.format(grmSD));
		}

		if (Math.abs(Effective_marker) > 0.0001) {
			Logger.printUserLog("Effective sample size is: " + df.format(Effective_sample));
		} else {
			Logger.printUserLog("Effective sample size is: " + dfE.format(Effective_sample));
		}
		
		if (Math.abs(Effective_marker) > 0.0001) {
			Logger.printUserLog("Effective number of genome segments is: " + df.format(Effective_marker));
		} else {
			Logger.printUserLog("Effective number of genome segments is: " + dfE.format(Effective_marker));
		}
		
		if(wgrmArgs.isGUI()) {
			PrintStream gui_file = null;
			gui_file = FileUtil.CreatePrintStream(wgrmArgs.getOutRoot()+".grm.gui");
			gui_file.println(df.format(Effective_sample) + "\t" + PI.size());
			gui_file.println(df.format(Effective_marker) + "\t" + allelefreq.length);
			gui_file.close();
		}

	}

	private double[] GRMScore(int idx1, int idx2) {
		double[] s = { 0, 0 };

		SubjectID subjectID1 = new SubjectID(pGM.getSample().get(idx1).getFamilyID(), pGM.getSample().get(idx1).getIndividualID());
		SubjectID subjectID2 = new SubjectID(pGM.getSample().get(idx2).getFamilyID(), pGM.getSample().get(idx2).getIndividualID());
		boolean idxFlag1 = subjectSet.contains(subjectID1);
		boolean idxFlag2 = subjectSet.contains(subjectID2);

		for (int i = 0; i < allelefreq.length; i++) {

			//inside control
			if (allelefreq[i][1] == Double.NaN || allelefreq[i][0] == 0 || allelefreq[i][1] == 0 || allelevar[i] == 0) continue;

			int g1 = pGM.getAdditiveScore(idx1, i);
			int g2 = pGM.getAdditiveScore(idx2, i);
			double m = allelefreq[i][1];
			if (g1 == Person.MissingGenotypeCode || g2 == Person.MissingGenotypeCode) {
				continue;
			} else {
				double de;
				if (wgrmArgs.isInbredList()) {
					if (idxFlag1 == true && idxFlag2 == true) {
						de = 4 * allelefreq[i][0] *  allelefreq[i][1];
					} if (idxFlag1 == false && idxFlag2 == false) {
						de = 2 * allelefreq[i][0] *  allelefreq[i][1];						
					} else {
						de = 2 * Math.sqrt(2) * allelefreq[i][0] *  allelefreq[i][1];						
					}
				} else {
					de = wgrmArgs.isInbred()? 4 * allelefreq[i][0] * allelefreq[i][1] : 2 * allelefreq[i][0] * allelefreq[i][1];
					if (wgrmArgs.isAdjVar()) de = allelevar[i];
				}
				s[0]++;
				s[1] += (g1 - 2 * m) * (g2 - 2 * m) / de;
			}
		}

		if (s[0] > 0) {
			s[1] /= s[0];
		} else {
			s[0] = 0;
			s[1] = 0;
		}
		return s;
	}

	private double[] GRMDomScore(int idx1, int idx2) {
		double[] s = { 0, 0 };

		for (int i = 0; i < allelefreq.length; i++) {
			//inside control
			if (allelefreq[i][1] == Double.NaN || allelefreq[i][0] == 0 || allelefreq[i][1] == 0 || allelevar[i] == 0) continue;

			double m = allelefreq[i][1];
			int g1 = pGM.getAdditiveScore(idx1, i);
			int g2 = pGM.getAdditiveScore(idx2, i);
			if (g1 == Person.MissingGenotypeCode || g2 == Person.MissingGenotypeCode) {
				continue;
			} else {
				s[0]++;

				double s1, s2;
				if (g1 == 0) {
					s1 = 0;
				} else if (g1 == 1) {
					s1 = 2 * m;
				} else {
					s1 = (4 * m - 2);
				}

				if (g2 == 0) {
					s2 = 0;
				} else if (g2 == 1) {
					s2 = 2 * m;
				} else {
					s2 = (4 * m - 2);
				}

				s[1] += (s1 - 2 * m * m) * (s2 - 2 * m * m) / (4 * m * m * (1 - m) * (1 - m));
			}
		}

		if (s[0] > 0) {
			s[1] /= s[0];
		} else {
			s[0] = 0;
			s[1] = 0;
		}

		return s;
	}
}
