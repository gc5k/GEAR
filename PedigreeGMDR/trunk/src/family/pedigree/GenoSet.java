package family.pedigree;

import java.util.Iterator;
import java.util.Set;
import java.util.TreeMap;

import publicAccess.PublicData;

public class GenoSet {

    TreeMap parentsGenoMap;
    TreeMap childrenGenoMap;
    TreeMap fullparentsGenoMap;
    TreeMap fullchildrenGenoMap;
    private int index;
    private int ungenotypedParents = 0;
    private int genotypedParents = 0;
    private int ungenotypedChildren = 0;
    private int genotypedChildren = 0;

    GenoSet(TreeMap parents, TreeMap children, int index) {
        fullparentsGenoMap = new TreeMap(parents);
        fullchildrenGenoMap = new TreeMap(children);
        parentsGenoMap = parents;
        childrenGenoMap = children;
        this.index = index;
        initial();
    }

    /**
     * kick the untyped individuals from the genotypeMap;
     */
    private void initial() {
        Set PG = parentsGenoMap.keySet();
        Iterator pit = PG.iterator();
        for (; pit.hasNext();) {
            genotypedParents += ((Integer) parentsGenoMap.get(pit.next())).intValue();
        }

        if (parentsGenoMap.containsKey(PublicData.MissingGenotype)) {
            ungenotypedParents = ((Integer) parentsGenoMap.get(PublicData.MissingGenotype)).intValue();
            parentsGenoMap.remove(PublicData.MissingGenotype);
        }
        genotypedParents -= ungenotypedParents;

        Set CG = childrenGenoMap.keySet();
        Iterator cit = CG.iterator();
        for (; cit.hasNext();) {
            genotypedChildren += ((Integer) childrenGenoMap.get(cit.next())).intValue();
        }

        if (childrenGenoMap.containsKey(PublicData.MissingGenotype)) {
            ungenotypedChildren = ((Integer) childrenGenoMap.get(PublicData.MissingGenotype)).intValue();
            childrenGenoMap.remove(PublicData.MissingGenotype);
        }
        genotypedChildren -= ungenotypedChildren;
    }

    public TreeMap getfullparentsGenoMap() {
        return fullparentsGenoMap;
    }

    /**
     * return typed parents genotypes'map
     * 
     * @return
     */
    public TreeMap getparentsGenoMap() {
        return parentsGenoMap;
    }

    public TreeMap getfullchildrenGenoMap() {
        return fullchildrenGenoMap;
    }

    /**
     * return typed children genotypes'map
     * 
     * @return
     */
    public TreeMap getchildrenGenoMap() {
        return (TreeMap) childrenGenoMap;
    }

    public int getNumParents() {
        return genotypedParents + ungenotypedParents;
    }

    public int getNumUntypedParents() {
        return ungenotypedParents;
    }

    public int getNumTypedParents() {
        return genotypedParents;
    }

    public int getNumChildren() {
        return genotypedChildren + ungenotypedChildren;
    }

    public int getNumUntypedChildren() {
        return ungenotypedChildren;
    }

    public int getNumTypedChildren() {
        return genotypedChildren;
    }

    public int getIndex() {
        return index;
    }
}