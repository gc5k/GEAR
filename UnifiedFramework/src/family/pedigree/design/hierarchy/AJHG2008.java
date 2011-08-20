package family.pedigree.design.hierarchy;

import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Random;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Map.Entry;

import org.apache.commons.lang3.ArrayUtils;

import util.NewIt;
import util.Sample;

import family.pedigree.genotype.GenoSet;

import family.RabinowitzLairdAlgorithm.AbstractGenoDistribution;
import family.RabinowitzLairdAlgorithm.lou.HeterozygousParent;
import family.RabinowitzLairdAlgorithm.lou.HomozygousParent;
import family.RabinowitzLairdAlgorithm.lou.ObservedParents;
import family.RabinowitzLairdAlgorithm.lou.UnobservedParents;
import family.mdr.MDRConstant;
import family.mdr.data.PersonIndex;
import family.pedigree.file.GMDRPhenoFile;
import family.pedigree.file.MapFile;
import family.pedigree.file.PedigreeFile;
import family.pedigree.genotype.BFamilyStruct;
import family.pedigree.genotype.BPerson;
import family.pedigree.phenotype.FamilyUnit;
import family.pedigree.phenotype.Subject;

public class AJHG2008 extends ChenBase {

	public AJHG2008(PedigreeFile ped, GMDRPhenoFile phe, MapFile map, long s, int pIdx, int[] cIdx, int m) {
		super(ped, phe, map, s, pIdx, cIdx, m);
		// TODO Auto-generated constructor stub
	}

	@Override
	protected void RevvingUp() {
		ChenBase.LineUpGenotypePhenotype lineup = new ChenBase.LineUpGenotypePhenotype();
		Hashtable<String, BFamilyStruct> Fam = PedData.getFamilyStruct();

		PersonTable.ensureCapacity(qualified_Sib);

		if (PhenoData != null)
			CovariateTable.ensureCapacity(qualified_Sib);

		// genotype = new byte[qualified_Sib][];
		status = new byte[qualified_Sib];

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
			FamilyUnit FamUnit = PhenoData == null ? null : PhenoData.getFamilyUnit(fi);
			String[] pi = fs.getPersonListSorted();
			int si = 0;
			ArrayList<BPerson> plist = NewIt.newArrayList();
			for (int i = 0; i < pi.length; i++) {
				if (!lineup.filter[c][i])
					continue;
				BPerson per = fs.getPerson(pi[i]);
				Subject sub = PhenoData == null ? null : FamUnit.getSubject(pi[i]);

				if (fs.hasAncestor(per)) {
					s_P.add(new PersonIndex(fs.getFamilyStructName(), pi[i], per));
					BPerson pseudoper = new BPerson(per);
					plist.add(pseudoper);
					status[s] = (byte) per.getAffectedStatus();
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
		AbstractGenoDistribution.rnd = new Random(seed);
		NonTransmitted(Fam);
	}

	public double[] getPermutedScore(boolean isNested) {

		if (isNested) {
			int c = 0;
			for (int i = 0; i < numSib.length; i++) {
				int[] si = Sample.SampleIndex(0, numSib[i] - 1, numSib[i]);
				for (int j = 0; j < si.length; j++) {
					PersonTable.get(c + j).setPermutedScore(score[c + si[j]]);
				}
				c += si.length;
			}
		} else {
			int[] idx = Sample.SampleIndex(0, PersonTable.size() - 1, PersonTable.size());
			for (int i = 0; i < idx.length; i++) {
				PersonTable.get(i).setPermutedScore(score[idx[i]]);
				// permuted_score[i] = score[idx[i]];
			}
		}

		return permuted_score;
	}

	private void NonTransmitted(Hashtable<String, BFamilyStruct> Fam) {

		Enumeration<String> perList;
		BFamilyStruct fam;
		BPerson per;
		BPerson pseudoper;
		String iid;
		String nontran_tran[];
		for (int i = 0; i < MapData.getMarkerNumber(); i++) {
			for (Entry<String, BFamilyStruct> entry : Fam.entrySet()) {
				fam = entry.getValue();
				perList = fam.getPersonList();
				GenoSet genoset = GenotypeSummary(fam, i);
				int numGenotypedParents = genoset.getNumTypedParents();
				AbstractGenoDistribution gDis;
				TreeSet<String> aSet = new TreeSet<String>();
				if (genoset.getNumParents() > 2) {
					System.err.println("Family " + entry.getKey() + " is not a nuclear family. It has more than 2 founders");
					System.exit(0);
				}
				if (numGenotypedParents == 2) {
					fam.countAllele(genoset.getchildrenGenoMap(), aSet);
					fam.countAllele(genoset.getparentsGenoMap(), aSet);
					if (aSet.size() > 4) {
						System.err.println("Marker " + i + " has more than 4 alleles.");
						System.exit(0);
					}
					gDis = new ObservedParents(genoset.getchildrenGenoMap(), genoset.getparentsGenoMap());
				} else {
					fam.countAllele(genoset.getchildrenGenoMap(), aSet);
					fam.countAllele(genoset.getparentsGenoMap(), aSet);
					if (numGenotypedParents == 1) {
						String PG = genoset.getparentsGenoMap().firstKey();
						if (!AbstractGenoDistribution.isHeterozygous(PG)) {// table
																			// 1
							if (aSet.size() > 3) {
								System.err.println("Marker " + i + " has more than 3 alleles with one parent is homozygous.");
								System.exit(0);
							}
							gDis = new HomozygousParent(genoset.getchildrenGenoMap(), genoset.getparentsGenoMap());
						} else {// table 2
							if (aSet.size() > 4) {
								System.err.println("Marker " + i + " has more than 4 alleles with one parent is heterozygous.");
								System.exit(0);
							}
							gDis = new HeterozygousParent(genoset.getchildrenGenoMap(), genoset.getparentsGenoMap());
						}
					} else {// table 3
						gDis = new UnobservedParents(genoset.getchildrenGenoMap());
					}
				}
				while (perList.hasMoreElements()) {
					iid = (String) perList.nextElement();
					per = fam.getPerson(iid);
					pseudoper = fam.getPseudoPerson(iid);
					if (fam.hasAncestor(per.getPersonID())) {
						continue;
					}
					nontran_tran = fam.getNonTransmitted(per.getGenotypeString(i), gDis);
					boolean f = nontran_tran[0].compareTo(MDRConstant.missingGenotype) == 0 ? false : true;
					if (f) {
						int a1 = Integer.parseInt(nontran_tran[0].substring(0, 1));
						int a2 = Integer.parseInt(nontran_tran[0].substring(1, 2));
						pseudoper.addMarker(f, a1, a2, i);
					} else {
						pseudoper.addMarker(f, 0, 0, i);
					}

				}
			}
		}
	}

	protected GenoSet GenotypeSummary(BFamilyStruct fs, int idx) {

		Boolean b = new Boolean(true);
		TreeMap<String, Integer> Ps;
		TreeMap<String, Integer> Ks;
		GenoSet gSet;

		Ps = NewIt.newTreeMap();
		Ks = NewIt.newTreeMap();
		Enumeration<String> perList = fs.getPersonList();
		while (perList.hasMoreElements()) {
			BPerson per = fs.getPerson((String) perList.nextElement());
			String genotype = per.getGenotypeString(idx);
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

}
