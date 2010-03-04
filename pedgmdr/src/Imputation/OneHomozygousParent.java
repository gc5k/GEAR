
package Imputation;

import java.util.Iterator;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

/**
 *
 * @author Guo-Bo Chen
 */
public class OneHomozygousParent extends AbstractImputation {
    TreeMap parentGenoMap;
    String parentgeno1;
    public OneHomozygousParent(TreeMap children, TreeMap parent) {
        super(children);
        parentGenoMap = new TreeMap(parent);
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
