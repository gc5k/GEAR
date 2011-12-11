package family.pedigree.design.hierarchy;

import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Random;

import admixture.parameter.Parameter;

import score.LinearRegression;
import score.LogisticRegression;
import test.Test;

import family.mdr.data.PersonIndex;
import family.mdr.partition.Partition;
import family.pedigree.file.MapFile;
import family.pedigree.file.PedigreeFile;
import family.pedigree.file.PhenotypeFile;
import family.pedigree.genotype.BFamilyStruct;
import family.pedigree.genotype.BPerson;
import family.pedigree.phenotype.FamilyUnit;
import family.pedigree.phenotype.Subject;

import util.NewIt;
import util.Sample;

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
	protected double[] status;
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
			boolean IsPhenotypeBinary = PedData.IsSixthColBinary();
			num_qualified = new int[Fam.size()][2];
			filter = new boolean[Fam.size()][];
			int c = 0;
			for (String fi : PedData.getFamListSorted()) {
				BFamilyStruct fs = Fam.get(fi);
				String[] pi = fs.getPersonListSorted();
				filter[c] = new boolean[pi.length];

				// filter_family

				if (PhenoData == null) {
					int cc = 0;
					for (int i = 0; i < pi.length; i++) {
						BPerson per = fs.getPerson(pi[i]);
						boolean hf = hardFilter(per);
						if (!hf) {
							filter[c][cc++] = hf;
							continue;
						}
						boolean f = true;
						if (IsPhenotypeBinary) {
							String s = per.getAffectedStatus();
							if (Parameter.status_shiftFlag) {
								f = (s.compareTo(Parameter.missing_phenotype) == 0) ? false
										: true;
							} else {
								f = (s.compareTo(Parameter.missing_phenotype) == 0 || s
										.compareTo("0") == 0) ? false : true;
							}
						} else {
							String s = per.getAffectedStatus();
							f = (s.compareTo(Parameter.missing_phenotype) == 0) ? true
									: false;
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
						boolean hf = hardFilter(per);
						if (!hf) {
							filter[c][cc++] = hf;
							continue;
						}
						boolean f = FamUnit.containsSubject(pi[i]);
						if (f) {
							Subject sub = FamUnit.getSubject(pi[i]);
							f = keep(fi, pi[i], per.getAffectedStatus(),
									sub.getTraits());
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

		protected boolean keep(String fid, String pid, String s,
				ArrayList<String> trait) {
			boolean f = true;
			boolean IsPhenotypeBinary = PedData.IsSixthColBinary();
			if (pheIdx == -1) {
				if (IsPhenotypeBinary) {
					if (Parameter.status_shiftFlag) {
						f = (s.compareTo(Parameter.missing_phenotype) == 0) ? false
								: true;
					} else {
						f = (s.compareTo(Parameter.missing_phenotype) == 0 || s
								.compareTo("0") == 0) ? false : true;
					}
				} else {
					f = (s.compareTo(Parameter.missing_phenotype) == 0) ? true
							: false;
				}
			} else {
				f = (trait.get(pheIdx).compareTo(Parameter.missing_phenotype) == 0) ? false
						: true;
			}

			if (covIdx == null) {
				return f;
			} else {
				for (int j = 0; j < covIdx.length; j++) {
					f = (trait.get(j).compareTo(Parameter.missing_phenotype) == 0) ? false
							: true;
				}
				return f;
			}
			/*
			 * further filtering if other conditions are applied.
			 */

		}

		protected boolean hardFilter(BPerson p) {
			boolean flag = true;

			if (Parameter.keepFlag) {
				flag = false;
				String fi = p.getFamilyID();
				String pi = p.getPersonID();
				for (int i = 0; i < Parameter.indKeep[0].length; i++) {
					if (fi.compareTo(Parameter.indKeep[0][i]) == 0
							&& pi.compareTo(Parameter.indKeep[1][i]) == 0) {
						flag = true;
						break;
					}
				}
			}
			if (Parameter.removeFlag) {
				String fi = p.getFamilyID();
				String pi = p.getPersonID();
				for (int i = 0; i < Parameter.ex_family[0].length; i++) {
					if (fi.compareTo(Parameter.ex_family[0][i]) == 0
							&& pi.compareTo(Parameter.ex_family[1][i]) == 0) {
						flag = false;
						break;
					}
				}
			} 
			if (flag && Parameter.keep_maleFlag) {
				return flag = p.getGender() == 1 ? true : false;
			}
			if (flag && Parameter.keep_femaleFlag) {
				return flag = p.getGender() == 2 ? true : false;
			}
			if (flag && Parameter.ex_nosexFlag) {
				return flag = (p.getGender() == 1 || p.getGender() == 2) ? true
						: false;
			}
			return flag;
		}
	}

	public ChenBase(PedigreeFile ped, PhenotypeFile phe, MapFile map, long s,
			int pIdx, int[] cIdx, int m) {
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
		if (PersonTable.size() < 10) {
			System.err.println("too few effective individuals ("
					+ PersonTable.size() + ") for the selected trait.");
			Test.LOG.append("too few effective individuals ("
					+ PersonTable.size() + ") for the selected trait.\n");
			Test.printLog();
			System.exit(0);
		}
		score = new double[PersonTable.size()];
		for (int i = 0; i < PersonTable.size(); i++) {
			if (pheIdx == -1) {
				if (PedData.IsSixthColBinary()) {
					score[i] = status[i] - Parameter.status_shift;
				} else {
					score[i] = status[i];
				}
			} else {

				if (PhenoData == null) {
					System.err.println("no phenotype file.");
					Test.LOG.append("no phenotype file.\n");
					Test.printLog();
					System.exit(0);
				}

				ArrayList<String> v = CovariateTable.get(i);
				String s = (String) v.get(pheIdx);
				score[i] = Double.parseDouble(s);
			}
		}
	}

	protected void buildScoreII() {
		if (PersonTable.size() < 10) {
			System.err.println("too few effective individuals ("
					+ PersonTable.size() + ") for the selected trait.");
			Test.LOG.append("too few effective individuals ("
					+ PersonTable.size() + ") for the selected trait.\n");
			Test.printLog();
			System.exit(0);
		}
		score = new double[PersonTable.size()];
		ArrayList<Double> T = NewIt.newArrayList();

		int method = 0;
		if(PedData.IsSixthColBinary()) {
			method = 1;
		}
		for (int i = 0; i < PersonTable.size(); i++) {
			double t = 0;
			if (PedData.IsSixthColBinary()) {
				t = status[i] - Parameter.status_shift;
			} else {
				t = status[i];
			}
			T.add(new Double(t));
		}

		double[][] X = null;

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
			LogisticRegression LogReg = new LogisticRegression(Y, X, true);
			LogReg.MLE();
			r = LogReg.getResiduals1();
		}

		System.arraycopy(r, 0, score, 0, r.length);
	}

	protected void buildScore(int pheIdx, int[] covIdx, int method) {
		if (PersonTable.size() < 10) {
			System.err.println("too few effective individuals ("
					+ PersonTable.size() + ") for the selected trait.");
			Test.LOG.append("too few effective individuals ("
					+ PersonTable.size() + ") for the selected trait.\n");
			Test.printLog();
			System.exit(0);
		}
		score = new double[PersonTable.size()];
		ArrayList<Double> T = NewIt.newArrayList();
		ArrayList<ArrayList<Double>> C = NewIt.newArrayList();

		for (int i = 0; i < PersonTable.size(); i++) {
			double t = 0;
			ArrayList<Double> c = NewIt.newArrayList();
			ArrayList<String> tempc = CovariateTable.get(i);
			if (pheIdx == -1) {
				if (PedData.IsSixthColBinary()) {
					t = status[i] - Parameter.status_shift;
				} else {
					t = status[i];
				}
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
	public double[] getStatus() {
		return status;
	}

	protected void generateScore() {

		if (PhenoData != null) {

			if(pheIdx == -1) {
				if (PedData.IsSixthColBinary()) {
					method = 1;
				} else {
					method = 0;
				}
			} else {

				for (int i = 0; i < PersonTable.size(); i++) {
					ArrayList<String> tempc = CovariateTable.get(i);
					String t = tempc.get(pheIdx);
					
					if(Parameter.status_shiftFlag) {
						if(t.compareTo("1") != 0 && t.compareTo("0")!= 0 && t.compareTo(Parameter.missing_phenotype) != 0) {
							method = 0;
							break;
						}
					} else {
						if(t.compareTo("2") != 0 && t.compareTo("1")!= 0 && t.compareTo("0") != 0 && t.compareTo(Parameter.missing_phenotype) != 0) {
							method = 0;
							break;
						}
					}
				}

			}

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
		ArrayList<Integer> g = null;
		if (Parameter.trgroupFlag) {
			g = Partition.TTPartition(PersonTable.size());
		} else if (Parameter.ttfileFlag) {
			g = Partition.ttFilePartition(PersonTable);
		} else if (Parameter.trsexFlag) {
			g = Partition.SexPartition(PersonTable);
		} else if (Parameter.cvFlag) {
			g = Partition.CVPartition(PersonTable.size());
		}

		for (int i = 0; i < PersonTable.size(); i++) {
			PersonIndex pi = PersonTable.get(i);
			pi.setGroup(g.get(i).intValue());
		}
	}

	@Override
	public void RecoverScore() {
		for (int i = 0; i < PersonTable.size(); i++) {
			PersonTable.get(i).setPermutedScore(score[i]);
		}
	}

	@Override
	public void recoverSeed() {
		rnd.setSeed(seed);
		Sample.setSeed(seed);
	}
}
