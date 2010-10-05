
package family.imputation;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import publicAccess.PublicData;
/**
 *
 * @author Guo-Bo Chen
 */
public class AbstractImputation {
    protected TreeMap childrenGenoMap;
    protected TreeSet alleleSet;
    protected TreeSet childrenalleleSet;
    protected TreeSet parentalleleSet;
    protected ArrayList parentGeno;

    public AbstractImputation(TreeMap children) {
        childrenGenoMap = new TreeMap(children);
        alleleSet = new TreeSet();
        childrenalleleSet = new TreeSet();
        parentalleleSet = new TreeSet();
        parentGeno = new ArrayList();
    }

    public String RandomAssign() {
        String geno = new String(PublicData.MissingGenotype);
        if (isParentGenotyped()) {
            String p1 = (String) parentGeno.get(0);
            String p2 = (String) parentGeno.get(1);
            char PG[][] = {{p1.charAt(0), p1.charAt(1)},
                {p2.charAt(0), p2.charAt(1)}
            };
            char allele[] = new char[2];
            int index1 = (Math.random() > 0.5) ? 0 : 1;
            int index2 = (Math.random() > 0.5) ? 0 : 1;
            if (PG[0][index1] <= PG[1][index2]) {
                allele[0] = PG[0][index1];
                allele[1] = PG[1][index2];
            } else {
                allele[1] = PG[0][index1];
                allele[0] = PG[1][index2];
            }
            geno = new String(allele);
        } else {
            if (childrenGenoMap.size() > 0) {
                Set GSet = childrenGenoMap.keySet();
                Iterator it = GSet.iterator();
                ArrayList Geno = new ArrayList();
                for (; it.hasNext();) {
                    Geno.add((String) it.next());
                }
                int index = (new Double(Math.random() * childrenGenoMap.size())).intValue();
                geno = new String((String) Geno.get(index));
            }
        }
        return geno;
    }

    protected boolean isParentGenotyped() {
        return parentGeno.size() == 2 ? true : false;
    }

    public void countChildrenAllele(Map Geno) {
        Set GSet = Geno.keySet();
        Iterator it = GSet.iterator();
        for (; it.hasNext();) {
            String g = (String) it.next();
            childrenalleleSet.add(g.substring(0, 1));
            childrenalleleSet.add(g.substring(1, 2));
        }
    }

    public void countParentAllele(TreeMap Geno) {
        Set GSet = Geno.keySet();
        Iterator it = GSet.iterator();
        for (; it.hasNext();) {
            String g = (String) it.next();
            parentalleleSet.add(g.substring(0, 1));
            parentalleleSet.add(g.substring(1, 2));
        }
    }

     /**
     * @param Geno
     *            Iterator the geno, and put alleles into alleleSet.
     */
    public void countAllele(Map Geno) {
        Set GSet = Geno.keySet();
        Iterator it = GSet.iterator();
        for (; it.hasNext();) {
            String g = (String) it.next();
            alleleSet.add(g.substring(0, 1));
            alleleSet.add(g.substring(1, 2));
        }
    }


    /**
     * @return the number of alleles within alleleSet.
     */
    public int numAllele() {
        return alleleSet.size();
    }

        /**
     * @param genotype
     * @return true if the genoytpe is heterozygous or false otherwise.
     */
    public static boolean isHeterozygous(String genotype) {
        return (genotype.charAt(0) != genotype.charAt(1)) ? true : false;
    }

    /**
     * @param genoMap
     * @return count the number of homozygous with in genoMap
     */
    public int numHomozygous(TreeMap genoMap) {
        int c = 0;
        Set gSet = genoMap.keySet();
        Iterator it = gSet.iterator();
        for (; it.hasNext();) {
            if (!isHeterozygous((String) it.next())) {
                c++;
            }
        }
        return c;
    }

    /**
     * @param genoMap
     * @return count the number of heterozygous within the genoMap
     */
    public int numHeterozygous(TreeMap genoMap) {
        int c = 0;
        Set gSet = genoMap.keySet();
        Iterator it = gSet.iterator();
        for (; it.hasNext();) {
            if (isHeterozygous((String) it.next())) {
                c++;
            }
        }
        return c;
    }


    /**
     * @return the whole set of alleles.
     */
    public TreeSet getAlleleSet() {
        return alleleSet;
    }
}
