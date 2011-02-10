package mdr.alternative;


import mdr.Combination;
import mdr.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import algorithm.CombinationGenerator;
import mdr.data.DataFile;

import mdr.Suite;
import publicAccess.PublicData;
import publicAccess.ToolKit;

/**
 *
 * @author Guo-Bo Chen
 */
public class GenotypeSearch extends AbstractSearch {

    public GenotypeSearch(DataFile dr, CombinationGenerator cg, int[] sI) {
        super(dr, cg, sI);
    }

    public void search(int or) {
        if (modelMap.containsKey(new Integer(or))) {
            return;
        }
        List sample = data.getSample();
        List com = (List) comGenerator.get(new Integer(or));
        HashMap submodelMap = new HashMap();
        for (Iterator e = com.iterator(); e.hasNext();) {
            HashMap subsets = new HashMap();
            String key = (String) e.next();
            int[] c = ToolKit.StringToIntArray(key);
            for (Iterator e1 = sample.iterator(); e1.hasNext();) {
                DataFile.Subject sub = (DataFile.Subject) e1.next();
                boolean genotyped = true;
                String genotype = new String();
                for (int j = 0; j < c.length; j++) {
                    String geno = (String) sub.getGenotype(c[j]);
                    if (geno.compareTo(PublicData.MissingGenotype) == 0) {
                        genotyped = false;
                        break;
                    }
                    if (genotype.length() == 0) {
                        genotype = genotype + geno;
                    } else {
                        genotype = genotype + PublicData.seperator + geno;
                    }
                }
                if (genotyped) {
                    ArrayList subset = (ArrayList) subsets.get(genotype);
                    if (subset == null) {
                        subset = new ArrayList();
                        subsets.put(genotype, subset);
                    }
                    subset.add(sub);
                }
            }

            Combination model = new Combination();
            Set genos = new TreeSet(subsets.keySet());
            for (Iterator e2 = genos.iterator(); e2.hasNext();) {
                String g = (String) e2.next();
                ArrayList subset = (ArrayList) subsets.get(g);
                Suite s = new Suite(subset, scrIdx.length);
                model.put(g, s);
            }
            submodelMap.put(key, model);
        }
        modelMap.put(new Integer(or), submodelMap);
    }
}
