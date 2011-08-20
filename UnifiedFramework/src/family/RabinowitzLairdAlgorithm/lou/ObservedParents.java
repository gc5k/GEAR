/**
 * 
 */
package family.RabinowitzLairdAlgorithm.lou;

import java.util.TreeMap;

import family.RabinowitzLairdAlgorithm.AbstractGenoDistribution;
import family.mdr.MDRConstant;

/**
 * Class for observing both parental genotypes.
 * 
 * @author Guobo Chen
 */
public class ObservedParents extends AbstractGenoDistribution {

	TreeMap<String, Integer> parentGenoMap;

	/**
	 * construct ObservedParents
	 * 
	 * @param children
	 *            genotypes of kids
	 * @param parents
	 *            genotypes of parents
	 */
	public ObservedParents(TreeMap<String, Integer> children, TreeMap<String, Integer> parents) {
		super(children);
		this.parentGenoMap = new TreeMap<String, Integer>(parents);
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
		char nontran[] = new char[2];
		String p1 = parentGeno.get(0);
		String p2 = parentGeno.get(1);
		char PG[][] = { { p1.charAt(0), p1.charAt(1) }, { p2.charAt(0), p2.charAt(1) } };
		char allele[] = new char[2];

		if (transmitted.compareTo(MDRConstant.missingGenotype) == 0) {// missing
																		// data
			int index1 = (rnd.nextFloat() > 0.5) ? 0 : 1;
			int index2 = (rnd.nextFloat() > 0.5) ? 0 : 1;
			if (PG[0][index1] <= PG[1][index2]) {
				allele[0] = PG[0][index1];
				allele[1] = PG[1][index2];
			} else {
				allele[1] = PG[0][index1];
				allele[0] = PG[1][index2];
			}
			char transallele[] = new char[2];
			if (PG[0][1 - index1] <= PG[1][1 - index2]) {
				transallele[0] = PG[0][1 - index1];
				transallele[1] = PG[1][1 - index2];
			} else {
				transallele[1] = PG[0][1 - index1];
				transallele[0] = PG[1][1 - index2];
			}
			String nontran_tran[] = { new String(allele), new String(transallele) };
			add(nontran_tran[0]);
			return nontran_tran;
		} else {
			allele[0] = transmitted.charAt(0);
			allele[1] = transmitted.charAt(1);
		}

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

		// for allele 1
		if (((L0[0][0] + L0[0][1]) > 0) && ((L0[1][0] + L0[1][1]) > 0)) {
			// found the allele in p1 and p2
			K0 = 2;
		} else if ((L0[0][0] + L0[0][1]) > 0) {// found the allele in p1;

			K0 = 0;
		} else {// found the allele in p2;

			K0 = 1;
		}

		// for allele 2
		if (((L1[0][0] + L1[0][1]) > 0) && ((L1[1][0] + L1[1][1]) > 0)) {
			// found the allele in p1 and p2
			K1 = 2;
		} else if ((L1[0][0] + L1[0][1]) > 0) {// found the allele in p1;

			K1 = 0;
		} else {// found the allele in p2;

			K1 = 1;
		}
		// System.out.println(C0 +" "+ C1 + " " + K0 + " " + K1);

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

		String nontran_tran[] = { new String(nontran), new String(transmitted) };
		add(nontran_tran[0]);

		return nontran_tran;
	}

	public void genotypeParents() {
		parentGeno.add(parentGenoMap.firstKey());
		parentGeno.add(parentGenoMap.lastKey());
	}
}
