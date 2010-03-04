package family.RabinowitzLairdAlgorithm;

import java.util.TreeMap;

import PublicAccess.PublicData;

/**
 * Class extends GenoDristribution Treat the situation of one heterozygous parent. However, it's still possible to
 * deduce the genotype of the other parent
 * 
 * @author Guobo Chen
 */
public class HomozygousParent extends AbstractGenoDistribution {

    TreeMap parentGenoMap;
    String parentgeno1;

    /**
     * @param child
     *            genotypes of kids
     * @param parent
     *            the genotype of the homozygous parent
     */
    public HomozygousParent(TreeMap child, TreeMap parent) {
        super(child);
        parentGenoMap = new TreeMap(parent);
        this.parentgeno1 = (String) parentGenoMap.firstKey();
        countAllele(childrenGenoMap);
        genotypeParents();
    }

    /**
     * Get nontransmitted genotypes If transmitted genotye is missing and can not be randomly assigned, the
     * nontransmitted genotype will be considered missing too.
     */
    public String[] getNontransmitted() {
        return null;
    }

    public String[] getNontransmitted(final String transmitted) {
        String nontran;
        String tran = new String(transmitted);
        String pgeno;

        if (transmitted.compareTo(PublicData.MissingGenotype) == 0) {
            tran = RandomAssign();
            if (tran.compareTo(PublicData.MissingGenotype) == 0) {
                String nontran_tran[] = {PublicData.MissingGenotype,
                    PublicData.MissingGenotype
                };
                return nontran_tran;
            }
        }
        if (childrenGenoMap.size() == 1) {// situation 1, 2

            pgeno = (String) childrenGenoMap.firstKey();
            nontran = new String(pgeno);
        } else {// situation 3, 4

            nontran = new String(tran.compareTo((String) childrenGenoMap.firstKey()) != 0 ? (String) childrenGenoMap.firstKey()
                    : (String) childrenGenoMap.lastKey());
        }
        add(nontran);
        String nontran_tran[] = {nontran, tran};
        return nontran_tran;
    }

    protected void genotypeParents() {
    }
}
