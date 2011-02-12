
package family.imputation;

import java.util.TreeMap;


/**
 *
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class GenotypedParents extends AbstractImputation{
    TreeMap<String, Integer> parentGenoMap;

    public GenotypedParents(TreeMap<String, Integer> children, TreeMap<String, Integer> parents) {
        super(children);
        this.parentGenoMap = new TreeMap<String, Integer>(parents);
        genotypeParents();
    }

    public void genotypeParents() {
        parentGeno.add((String) parentGenoMap.firstKey());
        parentGeno.add((String) parentGenoMap.lastKey());
    }
}
