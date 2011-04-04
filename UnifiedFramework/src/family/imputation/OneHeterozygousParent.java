
package family.imputation;

import java.util.TreeMap;
import java.util.TreeSet;

import util.NewIt;
/**
 *
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class OneHeterozygousParent extends AbstractImputation{
    TreeMap<String, Integer> parentGenoMap;
    String parentgeno1;
    public OneHeterozygousParent(TreeMap<String, Integer> children, TreeMap<String, Integer> parent) {
        super(children);
        parentGenoMap = new TreeMap<String, Integer>(parent);
        this.parentgeno1 = parentGenoMap.firstKey();
        countChildrenAllele(childrenGenoMap);
        countAllele(childrenGenoMap);
        countParentAllele(parentGenoMap);
        countAllele(parentGenoMap);
        genotypeParents();
    }

   /**
     * Genotype the other parent if possible.
     */
    protected void genotypeParents() {
        TreeSet<String> PG2 = NewIt.newTreeSet();
        if (numAllele() == 4) {// situation 16,17,18,19
            for (String it:alleleSet) {
                if (!parentgeno1.contains(it.substring(0, 1))) {
                    PG2.add(it);
                }
            }
        } else if (numAllele() == 3) {
            if (numHomozygous(childrenGenoMap) > 0) { // situation 8,9,10,11,12,13

                /**
                 * two step to genotype the second parent step 1: get the first allele from the homozygous genotype in
                 * children's geno set. step 2: get the second allele which is neither of the first parent's alleles.
                 */

                for (String it:childrenGenoMap.keySet()) {
                    if (!isHeterozygous(it)) {// step 1
                        PG2.add(it.substring(0, 1));
                        break;
                    }
                }

                for (String it:alleleSet) {
                    if (!parentgeno1.contains(it.substring(0, 1))) {// step 2

                        PG2.add(it);
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
            String parentgeno2 = new String(PG2.first() + PG2.last());
            parentGeno.add(parentgeno1);
            parentGeno.add(parentgeno2);
        }
    }
}
