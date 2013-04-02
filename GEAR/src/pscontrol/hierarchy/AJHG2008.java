package pscontrol.hierarchy;

import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Random;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Map.Entry;

import org.apache.commons.lang3.ArrayUtils;

import parameter.Parameter;

import family.pedigree.genotype.GenoSet;

import family.RabinowitzLairdAlgorithm.AbstractGenoDistribution;
import family.RabinowitzLairdAlgorithm.lou.HeterozygousParent;
import family.RabinowitzLairdAlgorithm.lou.HomozygousParent;
import family.RabinowitzLairdAlgorithm.lou.ObservedParents;
import family.RabinowitzLairdAlgorithm.lou.UnobservedParents;
import family.pedigree.Hukou;
import family.pedigree.PersonIndex;
import family.pedigree.file.MapFile;
import family.pedigree.file.PedigreeFile;
import family.pedigree.genotype.BFamilyStruct;
import family.pedigree.genotype.BPerson;
import gear.util.NewIt;

public class AJHG2008 extends ChenBase {

	private long seed = 2010;

	public AJHG2008(PedigreeFile ped, MapFile map) {
		super(ped, map);
	}
	
	public void RevvingUp(ArrayList<PersonIndex> pt) {
		Hashtable<String, BFamilyStruct> Fam = PedData.getFamilyStruct();
		Hashtable<String, BFamilyStruct> fam_has_sib = NewIt.newHashtable();
		PersonTable.ensureCapacity(qualified_Sib);

		ArrayList<PersonIndex> s_P = NewIt.newArrayList();

		ArrayList<Integer> SibIdx = NewIt.newArrayList();

		int c = 0;
		int s = 0;
		for (String fi : PedData.getFamListSorted()) {
			BFamilyStruct fs = Fam.get(fi);
			String[] pi = fs.getPersonListSorted();
			int si = 0;
			ArrayList<BPerson> plist = NewIt.newArrayList();
			for (int i = 0; i < pi.length; i++) {
				BPerson per = fs.getPerson(pi[i]);

				if (fs.hasAncestor(per)) {
					s_P.add(new PersonIndex(fs.getFamilyStructName(), pi[i],
							per, false, false));
					BPerson pseudoper = new BPerson(per);
					plist.add(pseudoper);
					s_P.add(new PersonIndex(fs.getFamilyStructName(), pseudoper
							.getPersonID(), pseudoper, true, false));
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

		numSib = ArrayUtils.toPrimitive(SibIdx.toArray(new Integer[0]));

		AbstractGenoDistribution.rnd = new Random(seed);
		NonTransmitted(fam_has_sib);
	}

	public void setSeed(long sd) {
		seed = sd;
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
					boolean f = g.compareTo(Parameter.INSTANCE.missingGenotype) != 0;
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

		TreeMap<String, Integer> Ps = NewIt.newTreeMap();
		TreeMap<String, Integer> Ks = NewIt.newTreeMap();
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

}
