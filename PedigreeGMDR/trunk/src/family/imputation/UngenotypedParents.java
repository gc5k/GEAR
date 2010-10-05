/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package family.imputation;

import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
/**
 *
 * @author Guo-Bo Chen
 */
public class UngenotypedParents extends AbstractImputation {
    public UngenotypedParents(TreeMap children) {
        super(children);
        countAllele(childrenGenoMap);
        genotypeParents();
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
}
