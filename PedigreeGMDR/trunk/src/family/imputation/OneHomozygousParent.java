
package family.imputation;

import java.util.TreeMap;

/**
 *
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class OneHomozygousParent extends AbstractImputation {
    TreeMap<String, Integer> parentGenoMap;
    String parentgeno1;
    public OneHomozygousParent(TreeMap<String, Integer> children, TreeMap<String, Integer> parent) {
        super(children);
        parentGenoMap = new TreeMap<String, Integer>(parent);
        this.parentgeno1 = (String) parentGenoMap.firstKey();
        countAllele(childrenGenoMap);
        genotypeParents();
    }
    /**
     * Genotype the other parent if possible.
     */
    protected void genotypeParents() {
    }
}
