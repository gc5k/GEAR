/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package family.imputation;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import util.NewIt;
/**
 *
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class UngenotypedParents extends AbstractImputation {
    public UngenotypedParents(TreeMap<String, Integer> children) {
        super(children);
        countAllele(childrenGenoMap);
        genotypeParents();
    }

    protected void genotypeParents() {
        ArrayList<TreeSet<String>> PG = NewIt.newArrayList();
        PG.add(new TreeSet<String>());
        PG.add(new TreeSet<String>());

        Set<String> GSet = childrenGenoMap.keySet();
        Iterator<String> it = GSet.iterator();
        String geno;
        if (numAllele() == 4 && childrenGenoMap.size() >= 3) {
            // three genotype is sufficient to recreat parents'genotypes
            geno = it.next();// the first genotype

            PG.get(0).add(new String(geno.substring(0, 1)));
            PG.get(1).add(new String(geno.substring(1, 2)));

            geno = it.next();
            if (!PG.get(0).contains(geno.substring(0, 1)) && !PG.get(0).contains(geno.substring(1, 2)) && !PG.get(1).contains(geno.substring(0, 1)) && !PG.get(1).contains(geno.substring(1, 2))) {// for situation 14
                // go to next genotype if this genotype does not match any allele being assigned in current stage.

                geno = it.next();
            }
            int index = 0;
            if (PG.get(0).contains(geno.substring(0, 1)) || PG.get(0).contains(geno.substring(1, 2))) {
                index = 1;
            }
            PG.get(index).add(new String((PG.get(1 - index).contains(geno.substring(0,
                    1))) ? geno.substring(1, 2) : geno.substring(0, 1)));
            geno = it.next();
            PG.get(1 - index).add(new String((PG.get(index).contains(geno.substring(0,
                    1))) ? geno.substring(1, 2) : geno.substring(0, 1)));
        }
        if (numAllele() == 3) {
            if (numHomozygous(childrenGenoMap) > 0) {
                for (; it.hasNext();) {
                    geno = (String) it.next();
                    if (!isHeterozygous(geno)) {
                        PG.get(0).add(new String(geno.substring(0, 1)));
                        PG.get(1).add(new String(geno.substring(0, 1)));
                    }
                }
                // System.out.println((String) PG.get(0).first());
                String HomoAllele = (String) PG.get(0).first();

                int index = 0;
                for (String al:getAlleleSet()) {
                    if (al.contains(HomoAllele.substring(0, 1))) {
                        continue;
                    }
                    PG.get(index++).add(new String(al));
                }
            }
        }
        if (numAllele() == 2) {
            if (numHomozygous(childrenGenoMap) == 2) {
                for (String al:getAlleleSet()) {
                    PG.get(0).add(new String(al));
                    PG.get(1).add(new String(al));
                }
            }
        }
        if (PG.get(0).size() == 2 && PG.get(1).size() == 2) {
            parentGeno.add(new String((String) PG.get(0).first() + (String) PG.get(0).last()));
            parentGeno.add(new String((String) PG.get(1).first() + (String) PG.get(1).last()));
        // System.out.println( "P1"+(String) parentGeno.firstElement() );
        // System.out.println( "P2"+(String) parentGeno.lastElement() );
        }
    }
}
