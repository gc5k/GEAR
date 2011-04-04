package mdr.alternative;

import mdr.*;
import java.util.ArrayList;
import java.util.HashMap;

import algorithm.CombinationGenerator;
import mdr.data.DataFile;

import mdr.Suite;
import publicAccess.PublicData;
import publicAccess.ToolKit;
import util.NewIt;

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
        ArrayList<String> com = comGenerator.get(new Integer(or));
        HashMap<String, Combination> submodelMap = NewIt.newHashMap();
        for (String key : com) {
            model = new Combination();
            SNPIndex = ToolKit.StringToIntArray(key);
            String c = new String();
            mergeSearch((ArrayList<DataFile.Subject>) data.getSample(), c, 0);
            submodelMap.put(key, model);
        }
        modelMap.put(new Integer(or), submodelMap);
    }

    public void mergeSearch(ArrayList<DataFile.Subject> subjects, String combination, int idxMarker) {
        if (idxMarker < SNPIndex.length) {
            HashMap<String, ArrayList<DataFile.Subject>> subsets = NewIt.newHashMap();
            for (DataFile.Subject sub: subjects) {
                String m = (String) sub.getGenotype(SNPIndex[idxMarker]);
                if (m.compareTo(PublicData.MissingGenotype) == 0) {
                    continue;
                }
                ArrayList<DataFile.Subject> subset = subsets.get(m);
                if (subset == null) {
                    subset = NewIt.newArrayList();
                    subsets.put(m, subset);
                }
                subset.add(sub);
            }

            for (String key : subsets.keySet()) {
                String com = new String();
                if (combination.length() == 0) {
                    com = key;
                } else {
                    com = combination + PublicData.seperator + key;
                }
                mergeSearch((ArrayList<DataFile.Subject>) subsets.get(key), com, idxMarker + 1);
            }
        } // if we've processed all attributes
        else {
            run(combination, subjects);
        }
    }

    public void run(String com, ArrayList<DataFile.Subject> subsample) {
        Suite s = new Suite(subsample, scrIdx.length);
        model.put(com, s);
    }

}
