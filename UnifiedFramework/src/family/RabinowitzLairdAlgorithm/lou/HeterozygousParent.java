package family.RabinowitzLairdAlgorithm.lou;

import java.util.ArrayList;
import java.util.TreeMap;
import java.util.TreeSet;

import family.RabinowitzLairdAlgorithm.AbstractGenoDistribution;
import family.mdr.arsenal.MDRConstant;

import util.NewIt;

/**
 * Class extends GenoDristribution Treat the situation of one heterozygous
 * parent. However, it's still possible to deduce the genotype of the other
 * parent
 * 
 * @author Guobo Chen, chenguobo@gmail.com
 */
public class HeterozygousParent extends AbstractGenoDistribution {

	TreeMap<String, Integer> parentGenoMap;
	String parentgeno1;

	/**
	 * @param children
	 *            genotypes of kids
	 * @param parent
	 *            the genotype of the heterozygous parent
	 */
	public HeterozygousParent(TreeMap<String, Integer> children, TreeMap<String, Integer> parent) {
		super(children);
		parentGenoMap = new TreeMap<String, Integer>(parent);
		this.parentgeno1 = parentGenoMap.firstKey();
		countChildrenAllele(childrenGenoMap);
		countAllele(childrenGenoMap);
		countParentAllele(parentGenoMap);
		countAllele(parentGenoMap);
		genotypeParents();
	}

	public String[] getNontransmitted() {
		return null;
	}

	/**
	 * Get nontransmitted genotypes If transmitted genotye is missing and can
	 * not be randomly assigned, the nontransmitted genotype will be considered
	 * missing too.
	 */
	public String[] getNontransmitted(final String transmitted) {
		String nontran = new String(MDRConstant.missingGenotype);
		String tran = new String(transmitted);

		if (transmitted.compareTo(MDRConstant.missingGenotype) == 0) {
			tran = RandomAssign();
			if (tran.compareTo(MDRConstant.missingGenotype) == 0) {
				String nontran_tran[] = { MDRConstant.missingGenotype, MDRConstant.missingGenotype };
				return nontran_tran;
			}
		}
		if (isParentGenotyped()) {
			String nontran_tran[] = produceNontransmitted(tran);
			return nontran_tran;
		} else if (childrenGenoMap.size() == 1) {
			if (!isHeterozygous(tran)) {// situation 1

				nontran = new String(parentgeno1);
			} else {
				if (tran.compareTo(parentgeno1) == 0) {// situation 2

					nontran = new String(tran);
				} else {// situation 3

					nontran = ExtractUniqueAllele2Genotype(parentgeno1, tran);
				}
			}
		} else if (childrenGenoMap.size() == 2) {
			if (numChildrenAllele() == 2) {// situation 4

				nontran = new String((tran.compareTo(childrenGenoMap.firstKey()) != 0) ? childrenGenoMap.firstKey() : childrenGenoMap.lastKey());

				if (!isHeterozygous(transmitted)) {// 4-1
					nontran = new String((transmitted.compareTo(childrenGenoMap.firstKey()) != 0) ? childrenGenoMap.firstKey() : childrenGenoMap
							.lastKey());
				} else {// 4-2
					double vk1 = (childrenGenoMap.get(childrenGenoMap.firstKey())).doubleValue();
					double vk2 = (childrenGenoMap.get(childrenGenoMap.lastKey())).doubleValue();
					double rd = rnd.nextFloat() * (vk1 + vk2);
					if (rd < vk1) {
						nontran = new String((String) childrenGenoMap.firstKey());
					} else {
						nontran = new String((String) childrenGenoMap.lastKey());
					}
				}

			} else if (numChildrenAllele() == 3) {
				if (!childrenGenoMap.containsKey(parentgeno1)) {// situation 5

					nontran = tran.compareTo(childrenGenoMap.firstKey()) != 0 ? childrenGenoMap.firstKey() : childrenGenoMap.lastKey();
				} else {// situation 14

					if (tran.compareTo(parentgeno1) == 0) {
						// code
						double rd = rnd.nextFloat();
						if (rd < 0.5) {
							nontran = ExtractUniqueAllele2Genotype(childrenGenoMap.firstKey(), childrenGenoMap.lastKey());
						} else {
							nontran = tran.compareTo(childrenGenoMap.firstKey()) != 0 ? childrenGenoMap.firstKey() : childrenGenoMap.lastKey();
						}
					} else {
						nontran = tran.compareTo(childrenGenoMap.firstKey()) != 0 ? childrenGenoMap.firstKey() : childrenGenoMap.lastKey();
					}
				}
			}
		} else {// situation 15

			if (tran.compareTo(parentgeno1) == 0) {
				double rd = rnd.nextFloat();
				ArrayList<String> GVec = NewIt.newArrayList();
				for (String geno : childrenGenoMap.keySet()) {
					if (tran.compareTo(geno) != 0) {
						GVec.add(geno);
					}
				}
				if (rd < 0.5) {
					nontran = GVec.get(0);
				} else {
					nontran = GVec.get(1);
				}
			} else {
				nontran = parentgeno1;
			}
		}
		add(nontran);
		String nontran_tran[] = { nontran, tran };
		return nontran_tran;
	}

	/**
	 * Genotype the other parent if possible.
	 */
	protected void genotypeParents() {
		TreeSet<String> PG2 = NewIt.newTreeSet();
		if (numAllele() == 4) {// situation 16,17,18,19

			for (String Callele : alleleSet) {
				if (!parentgeno1.contains(Callele.substring(0, 1))) {
					PG2.add(new String(Callele));
				}
			}
		} else if (numAllele() == 3) {
			if (numHomozygous(childrenGenoMap) > 0) { // situation
				// 8,9,10,11,12,13

				/**
				 * two step to genotype the second parent step 1: get the first
				 * allele from the homozygous genotype in children's geno set.
				 * step 2: get the second allele which is neither of the first
				 * parent's alleles.
				 */

				for (String CG : childrenGenoMap.keySet()) {
					if (!isHeterozygous(CG)) {// step 1

						PG2.add(CG.substring(0, 1));
						break;
					}
				}

				for (String Callele : alleleSet) {
					if (!parentgeno1.contains(Callele.substring(0, 1))) {// step
						// 2

						PG2.add(Callele);
						break;
					}
				}
			}
		} else if (numAllele() == 2) {
			if (numHomozygous(childrenGenoMap) == 2) {// situation 6, 7

				PG2 = new TreeSet<String>(getAlleleSet());
			}
		}

		if (PG2.size() > 0) {
			String parentgeno2 = PG2.first() + PG2.last();
			parentGeno.add(parentgeno1);
			parentGeno.add(parentgeno2);
		}
	}

	/**
	 * When both parental genotypes are available (the second parental genotype
	 * was deduced from previous steps)
	 * 
	 * @param transmitted
	 * @return String[]={nontransmitted, transmitted};
	 */
	public String[] produceNontransmitted(String transmitted) {
		char nontran[] = new char[2];
		// System.out.println((String) parentGeno.get(1));
		String p1 = parentGeno.get(0);
		String p2 = parentGeno.get(1);
		char[][] PG = { { p1.charAt(0), p1.charAt(1) }, { p2.charAt(0), p2.charAt(1) } };
		char allele[] = new char[2];

		allele[0] = transmitted.charAt(0);
		allele[1] = transmitted.charAt(1);

		int L0[][] = new int[2][2];
		int L1[][] = new int[2][2];
		int K0 = 0;
		int K1 = 0;
		int C0 = 0;
		int C1 = 0;
		for (int i = 0; i < PG.length; i++) {
			for (int j = 0; j < PG[i].length; j++) {
				if (allele[0] == PG[i][j]) {
					L0[i][j] = 1;
					C0++;
				}
				if (allele[1] == PG[i][j]) {
					L1[i][j] = 1;
					C1++;
				}
			}
		}

		if (((L0[0][0] + L0[0][1]) > 0) && ((L0[1][0] + L0[1][1]) > 0)) {// found
			// the
			// allele
			// in
			// p1
			// and
			// p2

			K0 = 2;
		} else if ((L0[0][0] + L0[0][1]) > 0) {// found the allele in p1;

			K0 = 0;
		} else {// found the allele in p2;

			K0 = 1;
		}

		if (((L1[0][0] + L1[0][1]) > 0) && ((L1[1][0] + L1[1][1]) > 0)) {// found
			// the
			// allele
			// in
			// p1
			// and
			// p2

			K1 = 2;
		} else if ((L1[0][0] + L1[0][1]) > 0) {// found the allele in p1;

			K1 = 0;
		} else {// found the allele in p2;

			K1 = 1;
		}
		/**
		 * All possible combination here are
		 * {(C0,C1)|(1,1),(1,2),(1,3),(2,2),(3,3)}
		 */
		if (C0 == 1 || C1 == 1) {// one allele meet one same allele in parent
			// genotype, the other more than one;

			if (C0 == 1) {
				K1 = 1 - K0;
				nontran[0] = (L0[K0][0] == 0) ? PG[K0][0] : PG[K0][1];
				nontran[1] = (L1[K1][0] == 0) ? PG[K1][0] : PG[K1][1];
			} else {
				K0 = 1 - K1;
				nontran[0] = (L0[K0][0] == 0) ? PG[K0][0] : PG[K0][1];
				nontran[1] = (L1[K1][0] == 0) ? PG[K1][0] : PG[K1][1];
			}
		} else {// c > 1 && c2 >1

			nontran[0] = (L0[0][0] == 0) ? PG[0][0] : PG[0][1];
			nontran[1] = (L1[1][0] == 0) ? PG[1][0] : PG[1][1];
		}

		if (nontran[1] < nontran[0]) {
			char temp = nontran[0];
			nontran[0] = nontran[1];
			nontran[1] = temp;
		}
		String tran_nontran[] = { new String(nontran), new String(transmitted) };
		return tran_nontran;
	}
}