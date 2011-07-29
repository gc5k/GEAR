package family.pedigree.design;

import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Map.Entry;

import util.NewIt;
import family.RabinowitzLairdAlgorithm.AbstractGenoDistribution;
import family.RabinowitzLairdAlgorithm.rabinowitz.Rabinowitz0;
import family.RabinowitzLairdAlgorithm.rabinowitz.Rabinowitz1;
import family.RabinowitzLairdAlgorithm.rabinowitz.Rabinowitz2;
import family.RabinowitzLairdAlgorithm.rabinowitz.Rabinowitz3;
import family.pedigree.genotype.FamilyStruct;
import family.pedigree.genotype.FamilyStructException;
import family.pedigree.genotype.GenoSet;
import family.pedigree.genotype.Person;

/**
*
* @author Guo-Bo Chen, chenguobo@gmail.com
*/
public class RLDriver {

	public RLDriver() {
		
	}

	public void TDT(Hashtable<String, FamilyStruct> familystructure, String[] markerInfor, int[] marker) {
//		Enumeration<String> fsList = familystructure.keys();
//		while (fsList.hasMoreElements()) {
		for(Entry<String, FamilyStruct> entry : familystructure.entrySet()) {
			String FID = entry.getKey();
			FamilyStruct fs = entry.getValue();
			TreeMap<String, Integer> Ps;
			TreeMap<String, Integer> Ks;
			GenoSet gSet = null;
			AbstractGenoDistribution gDis = null;
			for (int i = 0; i < marker.length; i++) {
				Ps = NewIt.newTreeMap();
				Ks = NewIt.newTreeMap();
				Enumeration<String> perList = fs.getPersonList();
				if (fs.getNumSibs() == 0 ) continue;
				while (perList.hasMoreElements()) {
					Person per = fs.getPerson(perList.nextElement());
					String genotype = per.getGenotype(marker[i]);
					if (fs.hasAncestor(per)) {
						if (Ks.containsKey(genotype)) {
							Integer c = ((Integer) Ks.get(genotype));
							Ks.put(genotype, ++c);
						} else {
							Ks.put(genotype, new Integer(1));
						}
					} else {
						if (Ps.containsKey(genotype)) {
							Integer c = ((Integer) Ps.get(genotype));
							Ps.put(genotype, ++c);
						} else {
							Ps.put(genotype, new Integer(1));
						}
					}
				}
				gSet = new GenoSet(Ps, Ks, marker[i]);

				TreeSet<String> aSet = NewIt.newTreeSet();
				if (gSet.getNumParents() > 2) {
					try {
						throw new FamilyStructException("Family " + FID
								+ " is not a nuclear family. It has more than 2 founders");
					} catch (FamilyStructException e) {
						e.printStackTrace();
					}
				}
				if (gSet.getNumTypedParents() == 2) {
					countAllele(gSet.getchildrenGenoMap(), aSet);
					countAllele(gSet.getparentsGenoMap(), aSet);
					gDis = new Rabinowitz0(gSet.getchildrenGenoMap(), gSet.getparentsGenoMap());
				} else {
					countAllele(gSet.getchildrenGenoMap(), aSet);
					countAllele(gSet.getparentsGenoMap(), aSet);
					if (gSet.getNumTypedParents() == 1) {
						String PG = (gSet.getparentsGenoMap()).firstKey();
						if (!AbstractGenoDistribution.isHeterozygous(PG)) {// table 1
							if (aSet.size() > 3) {
								try {
									throw new FamilyStructException("Marker " + markerInfor[marker[i]]
											+ " has more than 3 alleles with one parent is homozygous.");
								} catch (FamilyStructException e) {
									e.printStackTrace();
								}
							}
							gDis = new Rabinowitz1(gSet.getchildrenGenoMap(), gSet.getparentsGenoMap());
						} else {// table 2
							if (aSet.size() > 4) {
								try {
									throw new FamilyStructException("more than 4 alleles with one parent is heterozygous.");
								} catch (FamilyStructException e) {
									e.printStackTrace();
								}
							}
							gDis = new Rabinowitz2(gSet.getchildrenGenoMap(), gSet.getparentsGenoMap());
						}
					} else {// table 3
						gDis = new Rabinowitz3(gSet.getchildrenGenoMap());
					}
				}
				String[] controlGenotype = gDis.getNontransmitted();
				int index = 0;
				perList = fs.getPersonList();
				while (perList.hasMoreElements()) {
					String pi = perList.nextElement();
					Person per = (Person) fs.getPerson(pi);
					if (!fs.hasAncestor(per)) {
						continue;
					}
					per.setNonTransmittedGenotype(marker[i], controlGenotype[index++]);
				}
			}
		}
	}

	public void countAllele(TreeMap<String, Integer> Geno, Set<String> alleleSet) {

		for (String g : Geno.keySet()) {
			alleleSet.add(g.substring(0, 1));
			alleleSet.add(g.substring(1, 2));
		}
	}
}
