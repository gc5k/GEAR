package family.RabinowitzLairdAlgorithm.lou;

import java.util.TreeMap;

import family.RabinowitzLairdAlgorithm.AbstractGenoDistribution;
import family.mdr.arsenal.MDRConstant;


/**
 * Class extends GenoDristribution Treat the situation of one heterozygous
 * parent. However, it's still possible to deduce the genotype of the other
 * parent
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class HomozygousParent extends AbstractGenoDistribution {

	TreeMap<String, Integer> parentGenoMap;
	String parentgeno1;

	/**
	 * @param child
	 *            genotypes of kids
	 * @param parent
	 *            the genotype of the homozygous parent
	 */
	public HomozygousParent(TreeMap<String, Integer> child, TreeMap<String, Integer> parent) {
		super(child);
		parentGenoMap = new TreeMap<String, Integer>(parent);
		this.parentgeno1 = parentGenoMap.firstKey();
		countAllele(childrenGenoMap);
		genotypeParents();
	}

	/**
	 * Get nontransmitted genotypes If transmitted genotye is missing and can
	 * not be randomly assigned, the nontransmitted genotype will be considered
	 * missing too.
	 */
	public String[] getNontransmitted() {
		return null;
	}

	public String[] getNontransmitted(final String transmitted) {
		String nontran;
		String tran = new String(transmitted);
		String pgeno;

		if (transmitted.compareTo(MDRConstant.missingGenotype) == 0) {
			tran = RandomAssign();
			if (tran.compareTo(MDRConstant.missingGenotype) == 0) {
				String nontran_tran[] = { MDRConstant.missingGenotype, MDRConstant.missingGenotype };
				return nontran_tran;
			}
		}
		if (childrenGenoMap.size() == 1) {// situation 1, 2

			pgeno = (String) childrenGenoMap.firstKey();
			nontran = new String(pgeno);
		} else {// situation 3, 4

			nontran = tran.compareTo(childrenGenoMap.firstKey()) != 0 ? childrenGenoMap.firstKey() : childrenGenoMap.lastKey();
		}
		add(nontran);
		String nontran_tran[] = { nontran, tran };
		return nontran_tran;
	}

	protected void genotypeParents() {
		
	}
}
