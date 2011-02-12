
package family.RabinowitzLairdAlgorithm;

import java.util.ArrayList;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import publicAccess.PublicData;
import util.NewIt;

/**
 * An abstract class for Transmitted and Nontransmitted tables
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public abstract class AbstractGenoDistribution {

    protected TreeMap<String, Integer> childrenGenoMap;
    protected TreeMap<String, Integer> nontransmitted;
    protected TreeSet<String> alleleSet;
    protected TreeSet<String> childrenalleleSet;
    protected TreeSet<String> parentalleleSet;
    protected ArrayList<String> parentGeno;

    static public Random rnd;
    /**
     * Construct the GenoDistribution
     * 
     * @param child
     *            , is a TreeMap object containing numbers of different genotypes among the kids within the nuclear
     *            family.
     */
    public AbstractGenoDistribution(TreeMap<String, Integer> child) {
        childrenGenoMap = new TreeMap<String,Integer>(child);
        nontransmitted = NewIt.newTreeMap();
        alleleSet = NewIt.newTreeSet();
        childrenalleleSet = NewIt.newTreeSet();
        parentalleleSet = NewIt.newTreeSet();
        parentGeno = NewIt.newArrayList();
    }

    /**
     * It's a abstract method defined there
     * 
     * @param transmitted
     *            : the transmitted genotype of the person
     * @return String[]: nontrans_trans, the size is 2; nontran_tran[] = {nontrans, trans} trans will be the same
     *         reference of transmitted if transmitted is not missing and will be assigned a new referenc of transmitted
     *         genotype if it is missing.
     */
    public abstract String[] getNontransmitted(final String transmitted);

    public abstract String[] getNontransmitted();// RabinowitzProc

    /**
     * sometimes, the genotypes of parents are present; sometimes can be contructed from observed children genotype
     */
    protected abstract void genotypeParents();

    /**
     * @return true if both parental genotypes are available, or false otherwise.
     */
    protected boolean isParentGenotyped() {
        return parentGeno.size() == 2 ? true : false;
    }

    protected int getChildrenNum() {
        int sibs = 0;
        for (String it:childrenGenoMap.keySet()) {
            sibs += childrenGenoMap.get(it).intValue();
        }
        return sibs;
    }

    protected String[] getChildrenGenotypeArray() {
        Set<String> CSet = childrenGenoMap.keySet();
        String[] GenoSet = new String[CSet.size()];
        GenoSet = (String[]) CSet.toArray();
        return GenoSet;
    }

    protected void Produce(String[] control, TreeMap<String, Integer> cM, String[] genopool,
            double[] freq) {
        for (int i = 0; i < control.length; i++) {
            double rd = rnd.nextFloat();
            int index = 0;
            for (int j = 0; j < freq.length; j++) {
                index = j;
                if (rd <= freq[j]) {
                    break;
                }
            }
            control[i] = genopool[index];
            if (cM.containsKey(control[i])) {
                Integer c = ((Integer) cM.get(control[i]));
                c++;
                cM.put(control[i], c);
            } else {
                Integer c = new Integer(1);
                cM.put(new String(control[i]), c);
            }
        }
    }

    /**
     * @return a String with randomly assigned genotype. When parental genotypes are available, it's easy to randomly
     *         pick up an allele from each parent. If only sib genotypes are available, randomly pick up a genotype.
     */
    protected String RandomAssign() {
        String geno = new String(PublicData.MissingGenotype);
        if (isParentGenotyped()) {
            String p1 = (String) parentGeno.get(0);
            String p2 = (String) parentGeno.get(1);
            char PG[][] = {{p1.charAt(0), p1.charAt(1)},
                {p2.charAt(0), p2.charAt(1)}
            };
            char allele[] = new char[2];
            int index1 = (rnd.nextFloat() > 0.5) ? 0 : 1;
            int index2 = (rnd.nextFloat() > 0.5) ? 0 : 1;
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
                int index = (new Double(rnd.nextFloat() * childrenGenoMap.size())).intValue();
                geno = Geno.get(index);
            }
        }
        return geno;
    }

    public void countChildrenAllele(TreeMap<String, Integer> Geno) {

        for (String g:Geno.keySet()) {
            childrenalleleSet.add(g.substring(0, 1));
            childrenalleleSet.add(g.substring(1, 2));
        }
    }

    public TreeSet<String> getChildrenAlleleSet() {
        return childrenalleleSet;
    }

    public int numChildrenAllele() {
        return childrenalleleSet.size();
    }

    public void countParentAllele(TreeMap<String, Integer> Geno) {
        for (String g:Geno.keySet()) {
            parentalleleSet.add(g.substring(0, 1));
            parentalleleSet.add(g.substring(1, 2));
        }
    }

    public TreeSet<String> getParentAlleleSet() {
        return parentalleleSet;
    }

    public int numParentAllele() {
        return parentalleleSet.size();
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
     * @return the whole set of alleles.
     */
    public TreeSet<String> getAlleleSet() {
        return alleleSet;
    }

    /**
     * @return the number of alleles within alleleSet.
     */
    public int numAllele() {
        return alleleSet.size();
    }

    /**
     * Put a nontransmitted genotype into TreeMap nontransmitted and count the number of each produced nontransmitted
     * genotype.
     * 
     * @param genotype
     */
    public void add(String genotype) {
        if (nontransmitted.containsKey(genotype)) {
            Integer c = (Integer) nontransmitted.get(genotype);
            int v = c.intValue();
            v++;
            nontransmitted.put(genotype, new Integer(v));
        } else {
            nontransmitted.put(new String(genotype), new Integer(1));
        }
    }

    /**
     * @param genotype
     * @return return the count of the genotype.
     */
    public int get(String genotype) {
        if (nontransmitted.containsKey(genotype)) {
            return ((Integer) nontransmitted.get(genotype)).intValue();
        } else {
            return -1;
        }
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
     * genotype1 and genotype2 are different genotypes but with a same allele
     * 
     * @param genotype1
     * @param genotype2
     * @return a new genotype which is composed of the allele that are exclude. For exmaple, f(12, 13)=23.
     */
    public String ExtractUniqueAllele2Genotype(String genotype1,
            String genotype2) {
        String genotype;
        char allele1[] = {genotype1.charAt(0), genotype1.charAt(1)};
        char allele2[] = {genotype2.charAt(0), genotype2.charAt(1)};
        int index1 = 0, index2 = 0;
        if (allele1[0] == allele2[0]) {
            index1 = 1;
            index2 = 1;
        } else if (allele1[0] == allele2[1]) {
            index1 = 1;
            index2 = 0;
        } else if (allele1[1] == allele2[0]) {
            index1 = 0;
            index2 = 1;
        } else if (allele1[1] == allele2[1]) {
            index1 = 0;
            index2 = 0;
        }

        char allele[] = new char[2];
        if (genotype1.charAt(index1) <= genotype2.charAt(index2)) {
            allele[0] = genotype1.charAt(index1);
            allele[1] = genotype2.charAt(index2);
        } else {
            allele[0] = genotype2.charAt(index2);
            allele[1] = genotype1.charAt(index1);
        }

        genotype = new String(allele);
        return genotype;
    }
}
