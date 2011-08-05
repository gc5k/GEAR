package family.pedigree.design.hierarchy;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Random;

import admixture.parameter.Parameter;

import score.LinearRegression;
import score.LogisticRegression;

import family.pedigree.file.GMDRPhenoFile;
import family.pedigree.file.MapFile;
import family.pedigree.file.PedigreeFile;
import family.pedigree.file.SNP;
import family.pedigree.genotype.FamilyStruct;
import family.pedigree.genotype.Person;
import family.pedigree.phenotype.FamilyUnit;
import family.pedigree.phenotype.Subject;

import util.NewIt;

/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public abstract class ChenBase implements ChenInterface {

	protected long seed = 2011;
	public static Random rnd = new Random();
	protected MapFile MapData;
	protected PedigreeFile PedData;
	protected GMDRPhenoFile PhenoData;

	protected byte[][] genotype;
	protected int qualified_Unrelated;
	protected int qualified_Sib;
	protected int[] numSib;
	protected byte[] status;
	protected double[] score;
	protected double[] permuted_score;

	protected int pheIdx;
	protected int[] covIdx;
	protected int method;
	protected ArrayList<ArrayList<String>> CovariateTable;

	protected String[] scoreName = new String[1];
	protected ArrayList<PersonIndex> PersonTable;// The indexing file records
													// the

	protected int[] subsetMarker;

	public class PersonIndex {
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

	public class LineUpGenotypePhenotype {

		protected int[][] num_qualified;//
		protected boolean[][] filter;

		public LineUpGenotypePhenotype() {
			qualification();
		}

		private void qualification() {
			Hashtable<String, FamilyStruct> Fam = PedData.getFamilyStruct();
			num_qualified = new int[Fam.size()][2];
			filter = new boolean[Fam.size()][];
			int c = 0;
			for (String fi : PedData.getFamListSorted()) {
				FamilyStruct fs = Fam.get(fi);
				String[] pi = fs.getPersonListSorted();
				filter[c] = new boolean[pi.length];
				if (PhenoData != null && !PhenoData.containFamily(fi)) {  // if no phenotype for the whole family
					for (int i = 0; i < filter[c].length; i++) {
						filter[c][i] = false;
					}
				} else {
					FamilyUnit FamUnit = PhenoData.getFamilyUnit(fi);

					int cc = 0;
					for (int i = 0; i < pi.length; i++) {
						Person per = fs.getPerson(pi[i]);
						boolean f = FamUnit.containsSubject(pi[i]);
						if (f) {
							Subject sub = FamUnit.getSubject(pi[i]);
							f = filterItUp(fi, pi[i], (byte) per.getAffectedStatus(), sub.getTraits());
						}
						filter[c][cc++] = f;
						if (!f)
							continue;

						if (fs.hasAncestor(per)) {
							num_qualified[c][1]++;
						} else {
							num_qualified[c][0]++;
						}
					}
				}
				c++;
			}
			for (int i = 0; i < num_qualified.length; i++) {
				qualified_Unrelated += num_qualified[i][0];
				qualified_Sib += num_qualified[i][1];
			}
		}

		protected boolean filterItUp(String fid, String pid, byte s, ArrayList<String> trait) {
			boolean f = true;
			if (pheIdx == -1) {
				f = (s != (1 + Parameter.status_shift) && s != (2 + Parameter.status_shift)) ? false : true;
			} else {
				f = (trait.get(pheIdx).compareTo(Parameter.missing_phenotype) == 0) ? false : true;
			}
			if(PhenoData == null) {
				return f;
			}
			if (covIdx == null) {
				return f;
			} else {
				for (int j = 0; j < covIdx.length; j++) {
					f = (trait.get(j).compareTo(Parameter.missing_phenotype) == 0) ? false : true;
				}
				return f;
			}
			/*
			 * further filtering if other conditions are applied.
			 */

		}
	}

	public ChenBase(PedigreeFile ped, GMDRPhenoFile phe, MapFile map, long s, int pIdx, int[] cIdx, int m) {

		PedData = ped;
		PhenoData = phe;
		MapData = map;

		CovariateTable = NewIt.newArrayList();
		PersonTable = NewIt.newArrayList();
		pheIdx = pIdx;
		covIdx = cIdx;
		method = m;
		setSeed(s);
		RevvingUp();
		generateScore();
		
	}

	protected abstract void RevvingUp();

	private void fetchScore(int pheIdx) {
		double sum = 0;
		score = new double[PersonTable.size()];
		for (int i = 0; i < PersonTable.size(); i++) {
			if (pheIdx == -1) {
				score[i] = status[i] - 1;
			} else {
				try {
					if (PhenoData == null) {
						throw new Exception();
					} 
				} catch (Exception E) {
					System.err.println("no phenotype file");
				}
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

		for (int i = 0; i < PersonTable.size(); i++) {
			double t = 0;
			ArrayList<Double> c = NewIt.newArrayList();
			ArrayList<String> tempc = CovariateTable.get(i);
			if (pheIdx == -1) {// using affecting status as phenotype
				t = Parameter.status_shift == -1 ? status[i] : status[i] -1;
			} else {
				t = Double.parseDouble((String) tempc.get(pheIdx));
			}
			if (covIdx != null) {
				for (int j = 0; j < covIdx.length; j++) {
					c.add((Double.parseDouble((String) tempc.get(covIdx[j]))));
				}
				C.add(c);
			}
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

		if (CIndex != null) {
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
		ArrayList<SNP> snpList = MapData.getMarkerList();
		String[] m = new String[snpList.size()];
		for (int i = 0; i < snpList.size(); i++) {
			m[i] = snpList.get(i).getName();
		}
		return m;
	}

	@Override
	public double[] getPermutedScore(boolean nested) {
		return null;
	}

	@Override
	public double[][] getScore2() {
		double[][] s = new double[score.length][1];
		for (int i = 0; i < score.length; i++) {
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
		String[] mk = getMarkerName();
		for (int i = 0; i < mk.length; i++) {
			pedout.print(mk[i] + "\t");
		}
		pedout.println("status");

		for (int i = 0; i < genotype.length; i++) {
			for (int j = 0; j < genotype[i].length; j++) {
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

	private void generateScore() {
		if (method >= 0) {
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
	
	@Override
	public int getNumberMarker() {
		return MapData.getMarkerNumber();
	}

	@Override
	public MapFile getMapFile() {
		return MapData;
	}
}
