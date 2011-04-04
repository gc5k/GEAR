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
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class TreeSearch extends AbstractSearch {

    HashMap<String, Combination> models;

    public TreeSearch(DataFile d, CombinationGenerator cg, int[] sI) {
        super(d, cg, sI);
    }

    public void search(int or) {
        if(modelMap.containsKey(new Integer(or))) {
            return;
        }
        models = NewIt.newHashMap();
        order = or;
        String com = new String();
        int[] currIndex = new int[or];
        ArrayList<DataFile.Subject> subs = data.getSample();
        mergeSearch(0, subs, com, 1, currIndex);
        modelMap.put(new Integer(or), models);
    }

    public void mergeSearch(int fromIndex, ArrayList<DataFile.Subject> subjects, String combination, int depth, int[] currIndex) {
        if (depth <= order) {
            for (int or = fromIndex; or < (data.getMarkerNum() - order + depth); or++) {
                currIndex[depth - 1] = or;
                HashMap<String, ArrayList<DataFile.Subject>> genotypeSets = NewIt.newHashMap();
                for (DataFile.Subject r : subjects) {
                    String m = (String) r.getGenotype(or);
                    if (m.compareTo(PublicData.MissingGenotype) == 0) {
                        continue;
                    }
                    ArrayList<DataFile.Subject> subset = genotypeSets.get(m);
                    if (subset == null) {
                        subset = NewIt.newArrayList();
                        genotypeSets.put(m, subset);
                    }
                    subset.add(r);
                }

                for (String key:genotypeSets.keySet()) {
                    String com = new String();
                    if (combination.length() == 0) {
                        com = key;
                    } else {
                        com = combination + PublicData.seperator + key;
                    }
                    int[] cI = new int[currIndex.length];
                    System.arraycopy(currIndex, 0, cI, 0, currIndex.length);
                    mergeSearch((or + 1), (ArrayList<DataFile.Subject>) genotypeSets.get(key), com, (depth + 1), cI);
                }
            }
        } else {
            run(combination, subjects, currIndex);
        }
    }

    public void run(String com, ArrayList<DataFile.Subject> subsample, int[] cI) {
        Suite s = new Suite(subsample, scrIdx.length);
        String key = ToolKit.IntArrayToString(cI);
        Combination model = (Combination) models.get(key);
        if (model == null) {
            model = new Combination();
            models.put(key, model);
        }
        model.put(com, s);
    }
}
