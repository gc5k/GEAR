package family.pedigree.design.hierarchy;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
import java.util.Random;

import admixture.parameter.Parameter;

import score.LinearRegression;
import score.LogisticRegression;

import family.mdr.data.PersonIndex;
import family.pedigree.file.MapFile;
import family.pedigree.file.PedigreeFile;
import family.pedigree.file.PhenotypeFile;
import family.pedigree.genotype.BFamilyStruct;
import family.pedigree.genotype.BPerson;
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
	protected PhenotypeFile PhenoData;

	protected int qualified_Unrelated;
	protected int qualified_Sib;
	protected int[] numSib;
	protected byte[] status;
	protected double[] score;
	// protected double[] permuted_score;

	protected int pheIdx;
	protected int[] covIdx;
	protected int method;
	protected ArrayList<ArrayList<String>> CovariateTable;

	protected ArrayList<PersonIndex> PersonTable;// The indexing file records

	protected int[] subsetMarker;

	public class LineUpGenotypePhenotype {

		protected int[][] num_qualified;//
		protected boolean[][] filter;

		public LineUpGenotypePhenotype() {
			qualification();
		}

		private void qualification() {
			Hashtable<String, BFamilyStruct> Fam = PedData.getFamilyStruct();
			num_qualified = new int[Fam.size()][2];
			filter = new boolean[Fam.size()][];
			int c = 0;
			for (String fi : PedData.getFamListSorted()) {
				BFamilyStruct fs = Fam.get(fi);
				String[] pi = fs.getPersonListSorted();
				filter[c] = new boolean[pi.length];
				if (PhenoData == null) {
					int cc = 0;
					for (int i = 0; i < pi.length; i++) {
						BPerson per = fs.getPerson(pi[i]);
						int s = per.getAffectedStatus();
						boolean f = ((s + Parameter.status_shift) == 1 || (s + Parameter.status_shift) == 0) ? true : false;

						filter[c][cc++] = f;
						if (!f)
							continue;
						if (fs.hasAncestor(per)) {
							num_qualified[c][1]++;
						} else {
							num_qualified[c][0]++;
						}
					}
				} else if (!PhenoData.containFamily(fi)) {
					// if no phenotype for the whole family
					for (int i = 0; i < filter[c].length; i++) {
						filter[c][i] = false;
					}
				} else {
					FamilyUnit FamUnit = PhenoData.getFamilyUnit(fi);

					int cc = 0;
					for (int i = 0; i < pi.length; i++) {
						BPerson per = fs.getPerson(pi[i]);
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
				f = ((s + Parameter.status_shift) == 0 || (s + Parameter.status_shift) == 1) ? true : false;
			} else {
				f = (trait.get(pheIdx).compareTo(Parameter.missing_phenotype) == 0) ? false : true;
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

	public ChenBase(PedigreeFile ped, PhenotypeFile phe, MapFile map, long s, int pIdx, int[] cIdx, int m) {
		rnd.setSeed(Parameter.seed);
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

		group();

		CovariateTable = null;
	}

	protected abstract void RevvingUp();

	protected void fetchScore(int pheIdx) {
		double sum = 0;
		if(PersonTable.size() < 10) {
			System.err.println("too few effective individuals (" + PersonTable.size() + ") for the selected trait.");
			System.exit(0);
		}
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
	}

	protected void buildScoreII() {
		if(PersonTable.size() < 10) {
			System.err.println("too few effective individuals (" + PersonTable.size() + ") for the selected trait.");
			System.exit(0);
		}
		score = new double[PersonTable.size()];
		ArrayList<Double> T = NewIt.newArrayList();

		for (int i = 0; i < PersonTable.size(); i++) {
			double t = 0;
			t = status[i] + Parameter.status_shift;
			T.add(new Double(t));
		}

		double[][] X = null;

		double[][] Y = new double[T.size()][1];
		for (int i = 0; i < T.size(); i++) {
			Y[i][0] = ((Double) T.get(i)).doubleValue();
		}

		double[] r = null;

		LinearRegression LReg = new LinearRegression(Y, X, true);
		LReg.MLE();
		r = LReg.getResiduals1();

		System.arraycopy(r, 0, score, 0, r.length);
	}

	protected void buildScore(int pheIdx, int[] covIdx, int method) {
		if(PersonTable.size() < 10) {
			System.err.println("too few effective individuals (" + PersonTable.size() + ") for the selected trait.");
			System.exit(0);
		}
		score = new double[PersonTable.size()];
		ArrayList<Double> T = NewIt.newArrayList();
		ArrayList<ArrayList<Double>> C = NewIt.newArrayList();

		for (int i = 0; i < PersonTable.size(); i++) {
			double t = 0;
			ArrayList<Double> c = NewIt.newArrayList();
			ArrayList<String> tempc = CovariateTable.get(i);
			if (pheIdx == -1) {// using affecting status as phenotype
				t = status[i] + Parameter.status_shift;
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
	}

	@Override
	public void getPermutedScore(boolean nested) {

	}

	@Override
	public byte[] getStatus() {
		return status;
	}

	protected void generateScore() {
		if (PhenoData != null) {
			if (method >= 0) {
				buildScore(pheIdx, covIdx, method);
			} else {
				fetchScore(pheIdx);
			}
		} else {
			buildScoreII();
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

	public int SampleSize() {
		return PersonTable.size();
	}

	public ArrayList<PersonIndex> getSample() {
		return PersonTable;
	}

	private void group() {
		ArrayList<Integer> g = NewIt.newArrayList();
		for (int i = 0; i < PersonTable.size(); i++) {
			g.add(i % Parameter.cv);
		}
		Collections.shuffle(g, rnd);
		for (int i = 0; i < PersonTable.size(); i++) {
			PersonTable.get(i).setGroup(g.get(i));
		}
	}

	@Override
	public void RecoverScore() {
		for (int i = 0; i < PersonTable.size(); i++) {
			PersonTable.get(i).setPermutedScore(score[i]);
		}
	}
}
