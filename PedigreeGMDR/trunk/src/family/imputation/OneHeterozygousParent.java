
package family.imputation;

import java.util.Iterator;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
/**
 *
 * @author Guo-Bo Chen
 */
public class OneHeterozygousParent extends AbstractImputation{
    TreeMap parentGenoMap;
    String parentgeno1;
    public OneHeterozygousParent(TreeMap children, TreeMap parent) {
        super(children);
        parentGenoMap = new TreeMap(parent);
        this.parentgeno1 = (String) parentGenoMap.firstKey();
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
        TreeSet PG2 = new TreeSet();
        if (numAllele() == 4) {// situation 16,17,18,19

            Iterator it = alleleSet.iterator();
            for (; it.hasNext();) {
                String Callele = (String) it.next();
                if (!parentgeno1.contains(Callele.substring(0, 1))) {
                    PG2.add(new String(Callele));
                }
            }
        } else if (numAllele() == 3) {
            if (numHomozygous(childrenGenoMap) > 0) { // situation 8,9,10,11,12,13

                /**
                 * two step to genotype the second parent step 1: get the first allele from the homozygous genotype in
                 * children's geno set. step 2: get the second allele which is neither of the first parent's alleles.
                 */
                Set CSet = childrenGenoMap.keySet();
                Iterator it = CSet.iterator();
                for (; it.hasNext();) {
                    String CG = (String) it.next();
                    if (!isHeterozygous(CG)) {// step 1

                        PG2.add(CG.substring(0, 1));
                        break;
                    }
                }
                Iterator Ait = alleleSet.iterator();
                for (; Ait.hasNext();) {
                    String Callele = (String) Ait.next();
                    if (!parentgeno1.contains(Callele.substring(0, 1))) {// step 2

                        PG2.add(Callele);
                        break;
                    }
                }
            }
        } else if (numAllele() == 2) {
            if (numHomozygous(childrenGenoMap) == 2) {// situation 6, 7

                PG2 = new TreeSet(getAlleleSet());
            }
        }

        if (PG2.size() > 0) {
            String parentgeno2 = new String((String) PG2.first() + (String) PG2.last());
            parentGeno.add(parentgeno1);
            parentGeno.add(parentgeno2);
        }
    }
}
