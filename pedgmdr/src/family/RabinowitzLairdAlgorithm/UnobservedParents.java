package family.RabinowitzLairdAlgorithm;

import java.util.Iterator;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.ArrayList;

import PublicAccess.PublicData;

/**
 * Class extends GenoDistribution when neither parental genotype are available
 * 
 * @author Guobo Chen
 */
public class UnobservedParents extends AbstractGenoDistribution {

    /**
     * Construct UnobservedParents
     * 
     * @param child
     *            genotypes of observed kids
     */
    public UnobservedParents(TreeMap child) {
        super(child);
        countAllele(childrenGenoMap);
        genotypeParents();
    }

    public String[] getNontransmitted() {
        return null;
    }

    /**
     * Get nontransmitted genotypes If transmitted genotye is missing and can not be randomly assigned, the
     * nontransmitted genotype will be considered missing too.
     */
    public String[] getNontransmitted(final String transmitted) {
        String nontran = new String(PublicData.MissingGenotype);
        String tran = new String(transmitted);
        boolean flag = true;

        if (transmitted.compareTo(PublicData.MissingGenotype) == 0) {
            flag = false;
            tran = RandomAssign();
            if (tran.compareTo(PublicData.MissingGenotype) == 0) {
                String nontran_tran[] = {PublicData.MissingGenotype,
                    PublicData.MissingGenotype
                };
                return nontran_tran;
            }
        }
        if (isParentGenotyped()) {
            String nontran_tran[] = produceNontransmitted(tran);
            // System.out.println(nontran_tran[0]+"   "+nontran_tran[1]);
            return nontran_tran;
        } else if (childrenGenoMap.size() == 1) {// situation 1,2

            nontran = new String((String) childrenGenoMap.firstKey());
        } else if (childrenGenoMap.size() == 2) {
            if (numHomozygous(childrenGenoMap) > 0) {// situation 3

                nontran = (tran.compareTo((String) childrenGenoMap.firstKey()) != 0) ? (String) childrenGenoMap.firstKey() : (String) childrenGenoMap.lastKey();

                if (!isHeterozygous(transmitted)) {//3-1 
                    nontran = (transmitted.compareTo((String) childrenGenoMap.firstKey()) != 0) ? (String) childrenGenoMap.firstKey() : (String) childrenGenoMap.lastKey();
                } else {//3-2 
                    Set GSet = childrenGenoMap.keySet();
                    double vk1 = ((Integer) childrenGenoMap.get((String) childrenGenoMap.firstKey())).doubleValue();
                    double vk2 = ((Integer) childrenGenoMap.get((String) childrenGenoMap.lastKey())).doubleValue();
                    double rd = Math.random() * ( vk1 + vk2 );
                    if (rd < vk1) {
                        nontran = new String((String) childrenGenoMap.firstKey());
                    } else {
                        nontran = new String((String) childrenGenoMap.lastKey());
                    }
                }
            } else {// situation 6, 7

                nontran = new String((tran.compareTo((String) childrenGenoMap.firstKey()) != 0) ? (String) childrenGenoMap.firstKey() : (String) childrenGenoMap.lastKey());
            }
        } else {// situation 12

            Set GSet = childrenGenoMap.keySet();
            Iterator it = GSet.iterator();
            ArrayList Geno = new ArrayList();
            for (; it.hasNext();) {
                Geno.add((String) it.next());
            }
            int index = (new Double(Math.random() * childrenGenoMap.size())).intValue();
            nontran = (String) Geno.get(index);
        }
        add(nontran);

        String nontran_tran[] = {nontran, tran};
        // System.out.println(nontran_tran[0]+"   "+nontran_tran[1]);
        return nontran_tran;
    }

    protected void genotypeParents() {
        TreeSet PG[] = {new TreeSet(), new TreeSet()};

        Set GSet = childrenGenoMap.keySet();
        Iterator it = GSet.iterator();
        String geno;
        if (numAllele() == 4 && childrenGenoMap.size() >= 3) {
            // three genotype is sufficient to recreat parents'genotypes
            geno = (String) it.next();// the first genotype

            PG[0].add(new String(geno.substring(0, 1)));
            PG[1].add(new String(geno.substring(1, 2)));

            geno = (String) it.next();
            if (!PG[0].contains(geno.substring(0, 1)) && !PG[0].contains(geno.substring(1, 2)) && !PG[1].contains(geno.substring(0, 1)) && !PG[1].contains(geno.substring(1, 2))) {// for situation 14
                // go to next genotype if this genotype does not match any allele being assigned in current stage.

                geno = (String) it.next();
            }
            int index = 0;
            if (PG[0].contains(geno.substring(0, 1)) || PG[0].contains(geno.substring(1, 2))) {
                index = 1;
            }
            PG[index].add(new String((PG[1 - index].contains(geno.substring(0,
                    1))) ? geno.substring(1, 2) : geno.substring(0, 1)));
            geno = (String) it.next();
            PG[1 - index].add(new String((PG[index].contains(geno.substring(0,
                    1))) ? geno.substring(1, 2) : geno.substring(0, 1)));
        }
        if (numAllele() == 3) {
            if (numHomozygous(childrenGenoMap) > 0) {
                for (; it.hasNext();) {
                    geno = (String) it.next();
                    if (!isHeterozygous(geno)) {
                        PG[0].add(new String(geno.substring(0, 1)));
                        PG[1].add(new String(geno.substring(0, 1)));
                    }
                }
                // System.out.println((String) PG[0].first());
                String HomoAllele = (String) PG[0].first();
                TreeSet AS = getAlleleSet();
                Iterator iit = AS.iterator();
                int index = 0;
                for (; iit.hasNext();) {
                    String al = (String) iit.next();
                    if (al.contains(HomoAllele.substring(0, 1))) {
                        continue;
                    }
                    PG[index++].add(new String(al));
                }
            }
        }
        if (numAllele() == 2) {
            if (numHomozygous(childrenGenoMap) == 2) {
                TreeSet AS = getAlleleSet();
                Iterator iit = AS.iterator();
                for (; iit.hasNext();) {
                    String al = (String) iit.next();
                    PG[0].add(new String(al));
                    PG[1].add(new String(al));
                }
            }
        }
        if (PG[0].size() == 2 && PG[1].size() == 2) {
            parentGeno.add(new String((String) PG[0].first() + (String) PG[0].last()));
            parentGeno.add(new String((String) PG[1].first() + (String) PG[1].last()));
        // System.out.println( "P1"+(String) parentGeno.firstElement() );
        // System.out.println( "P2"+(String) parentGeno.lastElement() );
        }
    }

    public String[] produceNontransmitted(String transmitted) {
        char nontran[] = new char[2];
        String p1 = (String) parentGeno.get(0);
        String p2 = (String) parentGeno.get(1);
        char PG[][] = {{p1.charAt(0), p1.charAt(1)},
            {p2.charAt(0), p2.charAt(1)}
        };
        char allele[] = new char[2];
        boolean flag = true;

        allele[0] = transmitted.charAt(0);
        allele[1] = transmitted.charAt(1);

        int L0[][] = new int[2][2];
        int L1[][] = new int[2][2];
        int K0 = 0;
        int K1 = 0;
        int C0 = 0;
        int C1 = 0;
        for (int i = 0; i < PG.length; i++) {
            for (int j = 0; j < PG[i].length; j++) {
                if (allele[0] == PG[i][j]) {
                    L0[i][j] = 1;
                    C0++;
                }
                if (allele[1] == PG[i][j]) {
                    L1[i][j] = 1;
                    C1++;
                }
            }
        }

        if (((L0[0][0] + L0[0][1]) > 0) && ((L0[1][0] + L0[1][1]) > 0)) {// found the allele in p1 and p2

            K0 = 2;
        } else if ((L0[0][0] + L0[0][1]) > 0) {// found the allele in p1;

            K0 = 0;
        } else {// found the allele in p2;

            K0 = 1;
        }

        if (((L1[0][0] + L1[0][1]) > 0) && ((L1[1][0] + L1[1][1]) > 0)) {// found the allele in p1 and p2

            K1 = 2;
        } else if ((L1[0][0] + L1[0][1]) > 0) {// found the allele in p1;

            K1 = 0;
        } else {// found the allele in p2;

            K1 = 1;
        }

        /**
         * All possible combination here are {(C0,C1)|(1,1),(1,2),(1,3),(2,2),(3,3)}
         */
        if (C0 == 1 || C1 == 1) {// one allele meets one same allele in parent genotype, the other more than one;

            if (C0 == 1) {
                K1 = 1 - K0;
                nontran[0] = (L0[K0][0] == 0) ? PG[K0][0] : PG[K0][1];
                nontran[1] = (L1[K1][0] == 0) ? PG[K1][0] : PG[K1][1];
            } else {
                K0 = 1 - K1;
                nontran[0] = (L0[K0][0] == 0) ? PG[K0][0] : PG[K0][1];
                nontran[1] = (L1[K1][0] == 0) ? PG[K1][0] : PG[K1][1];
            }
        } else {// c > 1 && c2 >1

            nontran[0] = (L0[0][0] == 0) ? PG[0][0] : PG[0][1];
            nontran[1] = (L1[1][0] == 0) ? PG[1][0] : PG[1][1];
        }

        if (nontran[1] < nontran[0]) {
            char temp = nontran[0];
            nontran[0] = nontran[1];
            nontran[1] = temp;
        }

        String nontran_tran[] = {new String(nontran), new String(transmitted)};
        return nontran_tran;
    }
}