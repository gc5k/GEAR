
package family.imputation;

import java.util.ArrayList;
import java.util.TreeMap;
import java.util.TreeSet;

import publicAccess.PublicData;
import util.NewIt;
/**
 *
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class AbstractImputation {
    protected TreeMap<String, Integer> childrenGenoMap;
    protected TreeSet<String> alleleSet;
    protected TreeSet<String> childrenalleleSet;
    protected TreeSet<String> parentalleleSet;
    protected ArrayList<String> parentGeno;

    public AbstractImputation(TreeMap<String, Integer> children) {
        childrenGenoMap = new TreeMap<String, Integer>(children);
        alleleSet = NewIt.newTreeSet();
        childrenalleleSet = NewIt.newTreeSet();
        parentalleleSet = NewIt.newTreeSet();
        parentGeno = NewIt.newArrayList();
    }

    public String RandomAssign() {
        String geno = new String(PublicData.MissingGenotype);
        if (isParentGenotyped()) {
            String p1 = parentGeno.get(0);
            String p2 = parentGeno.get(1);
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
                ArrayList<String> Geno = NewIt.newArrayList();
                for (String it:childrenGenoMap.keySet()) {
                    Geno.add(it);
                }
                int index = (int) Math.random() * childrenGenoMap.size();
                geno = new String((String) Geno.get(index));
            }
        }
        return geno;
    }

    protected boolean isParentGenotyped() {
        return parentGeno.size() == 2 ? true : false;
    }

    public void countChildrenAllele(TreeMap<String, Integer> Geno) {
        for (String g:Geno.keySet()) {
            childrenalleleSet.add(g.substring(0, 1));
            childrenalleleSet.add(g.substring(1, 2));
        }
    }

    public void countParentAllele(TreeMap<String, Integer> Geno) {
        for (String g:Geno.keySet()) {
            parentalleleSet.add(g.substring(0, 1));
            parentalleleSet.add(g.substring(1, 2));
        }
    }

     /**
     * @param Geno
     *            Iterator the geno, and put alleles into alleleSet.
     */
    public void countAllele(TreeMap<String, Integer> Geno) {
        for (String g:Geno.keySet()) {
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
    public int numHomozygous(TreeMap<String, Integer> genoMap) {
        int c = 0;
        for (String g:genoMap.keySet()) {
            if (!isHeterozygous(g)) {
                c++;
            }
        }
        return c;
    }

    /**
     * @param genoMap
     * @return count the number of heterozygous within the genoMap
     */
    public int numHeterozygous(TreeMap<String, Integer> genoMap) {
        int c = 0;
        for (String g:genoMap.keySet()) {
            if (isHeterozygous(g)) {
                c++;
            }
        }
        return c;
    }


    /**
     * @return the whole set of alleles.
     */
    public TreeSet<String> getAlleleSet() {
        return alleleSet;
    }
}
