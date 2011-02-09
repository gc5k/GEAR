package mdrAlternative;

import mdr.Combination;
import mdr.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;

import algorithm.CombinationGenerator;
import mdr.data.DataFile;
import algorithm.Subdivision;

import mdr.Suite;
import publicAccess.PublicData;
import publicAccess.ToolKit;

/**
 *
 * @author Guo-Bo Chen
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
        models = new HashMap();
        order = or;
        String com = new String();
        int[] currIndex = new int[or];
        ArrayList subs = (ArrayList) data.getSample();
        mergeSearch(0, subs, com, 1, currIndex);
        modelMap.put(new Integer(or), models);
    }

    public void mergeSearch(int fromIndex, ArrayList subjects, String combination, int depth, int[] currIndex) {
        if (depth <= order) {
            for (int or = fromIndex; or < (data.getMarkerNum() - order + depth); or++) {
                currIndex[depth - 1] = or;
                HashMap genotypeSets = new HashMap();
                for (Iterator i = subjects.iterator(); i.hasNext();) {
                    DataFile.Subject r = (DataFile.Subject) i.next();
                    String m = (String) r.getGenotype(or);
                    if (m.compareTo(PublicData.MissingGenotype) == 0) {
                        continue;
                    }
                    ArrayList subset = (ArrayList) genotypeSets.get(m);
                    if (subset == null) {
                        subset = new ArrayList();
                        genotypeSets.put(m, subset);
                    }
                    subset.add(r);
                }
                Set Keys = genotypeSets.keySet();
                for (Iterator i = Keys.iterator(); i.hasNext();) {
                    String key = (String) i.next();
                    String com = new String();
                    if (combination.length() == 0) {
                        com = key;
                    } else {
                        com = combination + PublicData.seperator + key;
                    }
                    int[] cI = new int[currIndex.length];
                    System.arraycopy(currIndex, 0, cI, 0, currIndex.length);
                    mergeSearch((or + 1), (ArrayList) genotypeSets.get(key), com, (depth + 1), cI);
                }
            }
        } else {
            run(combination, subjects, currIndex);
        }
    }

    public void run(String com, ArrayList subsample, int[] cI) {
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
