package family.pedigree.design.hierarchy;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Random;

import score.LinearRegression;
import score.LogisticRegression;

import family.pedigree.file.GMDRPhenoFile;
import family.pedigree.file.GMDRPhenoFileException;
import family.pedigree.file.MDRPed;
import family.pedigree.file.MDRPedFileException;

import util.NewIt;

/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class ChenBase implements ChenInterface {
	
	protected long seed = 2011;
	public static Random rnd = new Random();
	protected MDRPed PedData;
	protected GMDRPhenoFile PhenoData;

	protected byte[][] genotype;
	protected int numUnrelated;
	protected int numSibs;
	protected int[] numSib;
	protected byte[] status;
	protected double[] score;
	protected double[] permuted_score;

	protected ArrayList<ArrayList<String>> CovariateTable;

	protected String[] scoreName = new String[1];
	protected ArrayList<PersonIndex> PersonTable;// The indexing file records the

	protected int[] subsetMarker;

	public static class PersonIndex {
		String FamilyID;
		String IndividualID;
		String Key;

		PersonIndex(String fid, String id) {
			FamilyID = new String(fid);
			IndividualID = new String(id);
			Key = new String(fid + "->" + id);
		}

		public String getFamilyID() {
			return FamilyID;
		}

		public String getIndividualID() {
			return IndividualID;
		}

		public String getKey() {
			return Key;
		}
	}

	public ChenBase(String ped, String phe) {

		PhenoData = new GMDRPhenoFile();
		PedData = new MDRPed();
		CovariateTable = NewIt.newArrayList();
		PersonTable = NewIt.newArrayList();
		
		RevvingUp(ped, phe);
	}

	protected void RevvingUp(String ped, String phe) {};

	private void fetchScore(int pheIdx) {
		double sum = 0;
		score = new double[PersonTable.size()];
		for (int i = 0; i < PersonTable.size(); i++) {
			if (pheIdx == -1) {
				score[i] = status[i] - 1;
			} else {
				ArrayList<String> v = CovariateTable.get(i);
				String s = (String) v.get(pheIdx);
				score[i] = Double.parseDouble(s);
				sum += score[i];
			}
		}
		if (pheIdx == -1) {
			scoreName[0] = new String("status");
		} else {
			scoreName[0] = new String(PhenoData.getTraitAtI(pheIdx));
		}
	}

	private void buildScore(int pheIdx, int[] covIdx, int method) {
		score = new double[PersonTable.size()];
		ArrayList<Double> T = NewIt.newArrayList();
		ArrayList<ArrayList<Double>> C = NewIt.newArrayList();
		ArrayList<PersonIndex> P = NewIt.newArrayList();

		for (int i = 0; i < PersonTable.size(); i++) {
			double t = 0;
			ArrayList<Double> c = NewIt.newArrayList();
			ArrayList<String> tempc = CovariateTable.get(i);
			if (pheIdx == -1) {// using affecting status as phenotype
				t = status[i];
			} else {
				t = Double.parseDouble((String) tempc.get(pheIdx));
			}
			if (covIdx != null) {
				for (int j = 0; j < covIdx.length; j++) {
					c.add((Double.parseDouble((String) tempc.get(covIdx[j]))));
				}
				C.add(c);
			}
			P.add(PersonTable.get(i));
			T.add(new Double(t));
		}

		double[][] X = null;
		if (covIdx != null) {
			X = new double[C.size()][covIdx.length];
			for (int i = 0; i < C.size(); i++) {
				ArrayList<Double> c = C.get(i);
				for (int j = 0; j < c.size(); j++) {
					X[i][j] = ((Double) c.get(j)).doubleValue();
				}
			}
		}

		double[][] Y = new double[T.size()][1];
		for (int i = 0; i < T.size(); i++) {
			Y[i][0] = ((Double) T.get(i)).doubleValue();
		}

		double[] r = null;
		if (method == 0) {
			LinearRegression LReg = new LinearRegression(Y, X, true);
			LReg.MLE();
			r = LReg.getResiduals1();
		} else {
			LogisticRegression LogReg1 = new LogisticRegression(Y, X, true);
			LogReg1.MLE();
			r = LogReg1.getResiduals1();
		}

		System.arraycopy(r, 0, score, 0, r.length);
		nameScore(pheIdx, covIdx, method);
	}

	private void nameScore(int PIndex, int[] CIndex, int method) {
		StringBuilder ln = new StringBuilder(300);
		if (method == 1) {
			ln.append("linear(");
		} else {
			ln.append("logistic(");
		}
		if (PIndex == -1) {
			ln.append("status->");
		} else {
			ln.append(PhenoData.getTraitAtI(PIndex));
		}
		ln.append("->");
		
		if(CIndex != null) {
			for (int i = 0; i < CIndex.length; i++) {
				ln.append(PhenoData.getTraitAtI(CIndex[i]) + ",");
			}
		} else {
			ln.append(",");
		}
		scoreName[0] = ln.toString();
	}

	@Override
	public byte[][] getGenotype() {
		return genotype;
	}

	@Override
	public String[] getMarkerName() {
		return PedData.getMarkerInformation(subsetMarker).toArray(new String[0]);
	}

	@Override
	public double[] getPermutedScore(boolean nested) {
		return null;
	}

	@Override
	public double[][] getScore2() {
		double[][] s = new double[score.length][1];
		for(int i = 0; i < score.length; i++) {
			s[i][0] = score[i];
		}
		return s;
	}

	@Override
	public String[] getScoreName() {
		return scoreName;
	}

	@Override
	public byte[] getStatus() {
		return status;
	}

	@Override
	public void print2MDRFormat(String f) {
		PrintWriter pedout = null;
		try {
			pedout = new PrintWriter(new File(f));
		} catch (FileNotFoundException e) {
			e.printStackTrace(System.err);
		}
		ArrayList<String> mk = PedData.getMarkerInformation();
		for(int i = 0; i < mk.size(); i++) {
			pedout.print(mk.get(i) + "\t");
		}
		pedout.println("status");
		
		for(int i = 0; i < genotype.length; i++) {
			for(int j = 0; j < genotype[i].length; j++) {
				pedout.print(genotype[i][j] + "\t");
			}
			pedout.println(status[i]);
		}
		pedout.close();
	}

	@Override
	public double[] getScore() {
		return score;
	}

	/**
	 * Initialize basic implementation of the genotype file.
	 * 
	 * @param Ped
	 *            the name of the pedigree file
	 * @throws IOException
	 */
	public void ParsePedFile(String ped) {
		File PedFile = new File(ped);
		try {
			PedData.Initial(PedFile);
			PedData.parseLinkage();
		} catch (MDRPedFileException e) {
			System.err.println("PedFile initial exception.");
			e.printStackTrace(System.err);
		} catch (IOException e) {
			System.err.println("Pedgree file initialization exception.");
			e.printStackTrace(System.err);
		}
	}

	/**
	 * Initialize basic implementation of the phenotype file.
	 * 
	 * @param Pheno
	 *            the name of the pedigree file
	 * @throws IOException
	 */
	public void ParsePhenoFile(String pheno) {
		File PheFile = new File(pheno);
		try {
			PhenoData.Initial(PheFile);
			PhenoData.parsePhenotype();
		} catch (GMDRPhenoFileException e) {
			System.err.println("Pheno file initialization exception.");
			e.printStackTrace(System.err);
		} catch (IOException e) {
			System.err.println("Pheno file initialization exception.");
			e.printStackTrace(System.err);
		}
	}

	public void SetChosenMarker(int[] mi) {
		subsetMarker = new int[mi.length];
		for (int i = 0; i < mi.length; i++) {
			subsetMarker[i] = mi[i];
		}
	}

	@Override
	public void generateScore(int pheIdx, int[] covIdx, int method) {
		if(method >= 0) {
			buildScore(pheIdx, covIdx, method);
		} else {
			fetchScore(pheIdx);
		}
		
	}

	@Override
	public void setSeed(long s) {
		seed = s;
		rnd.setSeed(s);
	}
}
