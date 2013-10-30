package family.pedigree.design.hierarchy;

import java.util.ArrayList;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Random;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Map.Entry;

import org.apache.commons.lang3.ArrayUtils;

import admixture.parameter.Parameter;

import score.LinearRegression;
import score.LogisticRegression;
import test.Test;
import util.NewIt;
import util.Sample;

import family.pedigree.genotype.GenoSet;

import family.RabinowitzLairdAlgorithm.AbstractGenoDistribution;
import family.RabinowitzLairdAlgorithm.lou.HeterozygousParent;
import family.RabinowitzLairdAlgorithm.lou.HomozygousParent;
import family.RabinowitzLairdAlgorithm.lou.ObservedParents;
import family.RabinowitzLairdAlgorithm.lou.UnobservedParents;
import family.mdr.arsenal.MDRConstant;
import family.mdr.data.PersonIndex;
import family.pedigree.file.MapFile;
import family.pedigree.file.PedigreeFile;
import family.pedigree.file.PhenotypeFile;
import family.pedigree.genotype.BFamilyStruct;
import family.pedigree.genotype.BPerson;
import family.pedigree.phenotype.FamilyUnit;
import family.pedigree.phenotype.Subject;

public class AJHG2008 extends ChenBase {

	public AJHG2008(PedigreeFile ped, PhenotypeFile phe, MapFile map, long s,
			int pIdx, int[] cIdx, int m) {
		super(ped, phe, map, s, pIdx, cIdx, m);
	}

	@Override
	protected void RevvingUp() {
		ChenBase.LineUpGenotypePhenotype lineup = new ChenBase.LineUpGenotypePhenotype();
		Hashtable<String, BFamilyStruct> Fam = PedData.getFamilyStruct();
		Hashtable<String, BFamilyStruct> fam_has_sib = NewIt.newHashtable();
		PersonTable.ensureCapacity(qualified_Sib);

		if (PhenoData != null)
			CovariateTable.ensureCapacity(qualified_Sib);

		status = new double[qualified_Sib];

		ArrayList<PersonIndex> s_P = NewIt.newArrayList();
		ArrayList<ArrayList<String>> s_C = NewIt.newArrayList();

		ArrayList<Integer> SibIdx = NewIt.newArrayList();

		int c = 0;
		int s = 0;
		for (String fi : PedData.getFamListSorted()) {
			if (lineup.num_qualified[c][1] == 0) {
				c++;
				continue;
			}
			BFamilyStruct fs = Fam.get(fi);
			FamilyUnit FamUnit = PhenoData == null ? null : PhenoData
					.getFamilyUnit(fi);
			String[] pi = fs.getPersonListSorted();
			int si = 0;
			ArrayList<BPerson> plist = NewIt.newArrayList();
			for (int i = 0; i < pi.length; i++) {
				if (!lineup.filter[c][i])
					continue;
				BPerson per = fs.getPerson(pi[i]);
				Subject sub = PhenoData == null ? null : FamUnit
						.getSubject(pi[i]);

				if (fs.hasAncestor(per.getPersonID())) {
					s_P.add(new PersonIndex(fs.getFamilyStructName(), pi[i],
							per, false, false));
					BPerson pseudoper = new BPerson(per);
					plist.add(pseudoper);
					s_P.add(new PersonIndex(fs.getFamilyStructName(), pseudoper
							.getPersonID(), pseudoper, true, false));
					status[s] = Double.parseDouble(per.getAffectedStatus());
					if (sub != null)
						s_C.add(sub.getTraits());
					si++;
					s++;
				}
			}
			if (si != 0) {
				SibIdx.add(new Integer(si));
				for (int i = 0; i < plist.size(); i++) {
					fs.addPerson(plist.get(i));
				}
				fam_has_sib.put(fs.getFamilyStructName(), fs);
			}
			c++;
		}
		PersonTable.addAll(s_P);
		if (PhenoData != null)
			CovariateTable.addAll(s_C);

		numSib = ArrayUtils.toPrimitive(SibIdx.toArray(new Integer[0]));

		int[] m = new int[MapData.getMarkerNumber()];
		for (int i = 0; i < MapData.getMarkerNumber(); i++) {
			m[i] = i;
		}
		int aff = 0;
		for (int i = 0; i < status.length; i++) {
			if ((status[i] - Parameter.status_shift) > 0) {
				aff++;
			}
		}
		
		System.err.println(s + " siblings " + "(" + aff + " affected "+ (s-aff)+" unaffected) from " + numSib.length + " families.");
		Test.LOG.append(s + " siblings " + "(" + aff + " affected "+ (s-aff)+" unaffected) from " + numSib.length + " families.\n");

		AbstractGenoDistribution.rnd = new Random(seed);
		NonTransmitted(fam_has_sib);
	}

	public void getPermutedScore(boolean isNested) {

		if (isNested) {
			int c = 0;
			for (int i = 0; i < numSib.length; i++) {
				if (rnd.nextBoolean()) {
					int[] si = Sample.SampleIndex(0, numSib[i] - 1, numSib[i]);
					for (int j = 0; j < si.length; j++) {
						PersonTable.get(c + j * 2).setPermutedScore(
								score[c + j * 2 + 1]);
						PersonTable.get(c + j * 2 + 1).setPermutedScore(
								score[c + j * 2]);
					}
				}
				c += numSib[i] * 2;
			}
		} else {
			int[] idx = Sample.SampleIndex(0, PersonTable.size() - 1,
					PersonTable.size());
			for (int i = 0; i < idx.length; i++) {
				PersonTable.get(i).setPermutedScore(score[idx[i]]);
			}
		}

	}

	private void NonTransmitted(Hashtable<String, BFamilyStruct> Fam) {

		Enumeration<String> perList;
		BFamilyStruct fam;
		BPerson per;
		BPerson pseudoper;
		String iid;

		for (int i = 0; i < MapData.getMarkerNumber(); i++) {
			for (Entry<String, BFamilyStruct> entry : Fam.entrySet()) {
				fam = entry.getValue();
				perList = fam.getPersonList();
				GenoSet genoset = GenotypeSummary(fam, i);
				int numGenotypedParents = genoset.getNumTypedParents();
				AbstractGenoDistribution gDis;
				TreeSet<String> aSet = new TreeSet<String>();
				if (numGenotypedParents == 2) {
					fam.countAllele(genoset.getchildrenGenoMap(), aSet);
					fam.countAllele(genoset.getparentsGenoMap(), aSet);
					gDis = new ObservedParents(genoset.getchildrenGenoMap(),
							genoset.getparentsGenoMap());
				} else {
					fam.countAllele(genoset.getchildrenGenoMap(), aSet);
					fam.countAllele(genoset.getparentsGenoMap(), aSet);
					if (numGenotypedParents == 1) {
						String PG = genoset.getparentsGenoMap().firstKey();
						if (!AbstractGenoDistribution.isHeterozygous(PG)) {
							// table 1
							gDis = new HomozygousParent(
									genoset.getchildrenGenoMap(),
									genoset.getparentsGenoMap());
						} else {// table 2
							gDis = new HeterozygousParent(
									genoset.getchildrenGenoMap(),
									genoset.getparentsGenoMap());
						}
					} else {// table 3
						gDis = new UnobservedParents(
								genoset.getchildrenGenoMap());
					}
				}
				perList = fam.getPersonList();
				while (perList.hasMoreElements()) {
					iid = perList.nextElement();
					if (iid.contains("ajhg2008") || !fam.hasAncestor(iid)) {
						continue;
					}
					per = fam.getPerson(iid);
					StringBuffer sb = new StringBuffer(iid);
					sb.append("ajhg2008");
					pseudoper = fam.getPerson(sb.toString());
					String g = per.getBiAlleleGenotypeString(i);
					boolean f = g.compareTo(MDRConstant.missingGenotype) != 0;
					String[] nontran_tran = new String[2];
					if (f) {
						nontran_tran = fam.getNonTransmitted(g, gDis);
						int a1 = Integer.parseInt(nontran_tran[0].substring(0,
								1));
						int a2 = Integer.parseInt(nontran_tran[0].substring(1,
								2));
						pseudoper.addMarker(f, a1, a2, i);
					} else {
						pseudoper.addMarker(f, 0, 0, i);
					}
				}
			}
		}
	}

	protected GenoSet GenotypeSummary(BFamilyStruct fs, int idx) {

		TreeMap<String, Integer> Ps;
		TreeMap<String, Integer> Ks;

		Ps = NewIt.newTreeMap();
		Ks = NewIt.newTreeMap();
		Enumeration<String> perList = fs.getPersonList();
		while (perList.hasMoreElements()) {
			String iid = perList.nextElement();
			if (iid.contains("ajhg2008"))
				continue;
			BPerson per = fs.getPerson(iid);
			String genotype = per.getBiAlleleGenotypeString(idx);
			if (fs.hasAncestor(per.getPersonID())) {
				if (Ks.containsKey(genotype)) {
					Integer c = ((Integer) Ks.get(genotype));
					Ks.put(genotype, ++c);
				} else {
					Integer c = new Integer(1);
					Ks.put(new String(genotype), c);
				}
			} else {
				if (Ps.containsKey(genotype)) {
					Integer c = ((Integer) Ps.get(genotype));
					Ps.put(genotype, ++c);
				} else {
					Integer c = new Integer(1);
					Ps.put(new String(genotype), c);
				}
			}
		}
		return new GenoSet(Ps, Ks, idx);
	}

	protected void fetchScore(int pheIdx) {
		score = new double[PersonTable.size()];
		for (int i = 0; i < PersonTable.size(); i += 2) {
			if (pheIdx == -1) {
				if (PedData.IsSixthColBinary()) {
					score[i] = status[i / 2] - Parameter.status_shift;
				} else {
					score[i] = status[i / 2];
				}
				score[i + 1] = -1 * score[i];
			} else {
				if (PhenoData == null) {
					System.err.println("no phenotype file.");
					Test.LOG.append("no phenotype file.\n");
					Test.printLog();
					System.exit(0);
				}
				ArrayList<String> v = CovariateTable.get(i / 2);
				String s = v.get(pheIdx);
				score[i] = Double.parseDouble(s);
				score[i + 1] = -1 * score[i];
			}
		}
	}

	@Override
	protected void buildScore(int pheIdx, int[] covIdx, int method) {
		score = new double[PersonTable.size()];
		ArrayList<Double> T = NewIt.newArrayList();
		HashMap<String, Integer> PheCat = NewIt.newHashMap();
		ArrayList<ArrayList<Double>> C = NewIt.newArrayList();

		for (int i = 0; i < PersonTable.size(); i += 2) {
			double t = 0;
			String phe;
			ArrayList<Double> c = NewIt.newArrayList();
			ArrayList<String> tempc = CovariateTable.get(i / 2);
			if (pheIdx == -1) {// using affecting status as phenotype
				t = status[i / 2] - Parameter.status_shift;
			} else {
				t = Double.parseDouble((String) tempc.get(pheIdx));
			}
			phe = Double.toString(t);
			if(PheCat.containsKey(phe)) {
				PheCat.put(phe, 1);
			} else {
				PheCat.put(phe, 1);
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
		if (PheCat.size() > 1) {
			LinearRegression LReg = new LinearRegression(Y, X, true);
			LReg.MLE();
			r = LReg.getResiduals1();			
			for (int i = 0; i < r.length; i++) {
				score[i * 2] = r[i];
				score[i * 2 + 1] = -1 * r[i];
			}

		} else {
			for (int i = 0; i < r.length; i++) {
				score[i * 2] = Y[i][0];
				score[i * 2 + 1] = -1 * Y[i][0];
			}
		}
	}

	@Override
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

		for (int i = 0; i < PersonTable.size() / 2; i++) {
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
			score[i * 2] = Y[i][0];
			score[i * 2 + 1] = -1 * Y[i][0];
		}

//		double[] r = null;
//		if (method == 0) {
//			LinearRegression LReg = new LinearRegression(Y, X, true);
//			LReg.MLE();
//			r = LReg.getResiduals1();
//		} else {
//			LogisticRegression LogReg1 = new LogisticRegression(Y, X, true);
//			LogReg1.MLE();
//			r = LogReg1.getResiduals1();
//		}
//
//		for (int i = 0; i < r.length; i++) {
//			score[i * 2] = r[i];
//			score[i * 2 + 1] = -1 * r[i];
//		}
	}
}
