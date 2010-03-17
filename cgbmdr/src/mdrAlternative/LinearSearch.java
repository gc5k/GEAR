package mdrAlternative;

import mdr.Combination;
import mdr.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import algorithm.CombinationGenerator;
import mdr.data.DataFile;

import mdr.Suite;
import publicAccess.PublicData;
import publicAccess.ToolKit;

/**
 *
 * @author Guo-Bo Chen
 */
public class LinearSearch extends AbstractSearch {
    
    private int[] SNPIndex;
    private Combination model;

    public LinearSearch(DataFile dr, CombinationGenerator cg, int[] sI) {
        super(dr, cg, sI);
    }

    public void search(int or) {
        if (modelMap.containsKey(new Integer(or))) {
            return;
        }
        List com = (List) comGenerator.get(or);
        HashMap submodelMap = new HashMap();
        for (Iterator e = com.iterator(); e.hasNext();) {
            model = new Combination();
            String key = (String) e.next();
            SNPIndex = ToolKit.StringToIntArray(key);
            String c = new String();
            mergeSearch((ArrayList) data.getSample(), c, 0);
            submodelMap.put(key, model);
        }
        modelMap.put(new Integer(or), submodelMap);
    }

    public void mergeSearch(ArrayList subjects, String combination, int idxMarker) {
        if (idxMarker < SNPIndex.length) {
            HashMap subsets = new HashMap();
            for (Iterator i = subjects.iterator(); i.hasNext();) {
                DataFile.Subject sub = (DataFile.Subject) i.next();
                String m = (String) sub.getGenotype(SNPIndex[idxMarker]);
                if (m.compareTo(PublicData.MissingGenotype) == 0) {
                    continue;
                }
                ArrayList subset = (ArrayList) subsets.get(m);
                if (subset == null) {
                    subset = new ArrayList();
                    subsets.put(m, subset);
                }
                subset.add(sub);
            }
            Set keys = subsets.keySet();
            for (Iterator i = keys.iterator(); i.hasNext();) {
                String key = (String) i.next();
                String com = new String();
                if (combination.length() == 0) {
                    com = key;
                } else {
                    com = combination + PublicData.seperator + key;
                }
                mergeSearch((ArrayList) subsets.get(key), com, idxMarker + 1);
            }
        } // if we've processed all attributes
        else {
            run(combination, subjects);
        }
    }

    public void run(String com, ArrayList subsample) {
        Suite s = new Suite(subsample, scrIdx.length);
        model.put(com, s);
    }

}
