
package family.imputation;

import java.util.Iterator;
import java.util.TreeMap;
import java.util.TreeSet;

import publicAccess.PublicData;

/**
 *
 * @author Guo-Bo Chen
 */
public class GenotypedParents extends AbstractImputation{
    TreeMap parentGenoMap;

    public GenotypedParents(TreeMap children, TreeMap parents) {
        super(children);
        this.parentGenoMap = new TreeMap(parents);
        genotypeParents();
    }

    public void genotypeParents() {
        parentGeno.add((String) parentGenoMap.firstKey());
        parentGeno.add((String) parentGenoMap.lastKey());
    }
}
