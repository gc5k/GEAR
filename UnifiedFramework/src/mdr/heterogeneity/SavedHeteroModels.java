package mdr.heterogeneity;

import mdr.heterogeneity.HeteroCombination;
import java.util.Iterator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

/**
 *
 * @author Guo-Bo Chen
 */
public class SavedHeteroModels {

    HashMap<String, HeteroCombination> savedModels = new HashMap();
    HashMap<String, Integer> modelCount = new HashMap();
    HashMap<Integer, String> modelKeys = new HashMap();//keep accordingly the key of the best model to each trait.

    public SavedHeteroModels() {}

    public void save(String key, HeteroCombination model) {
        if( !savedModels.containsKey(key)) {
            savedModels.put(key, model);
        }
        Integer value = (Integer) modelCount.get(key);
        if( value == null ) {
            modelCount.put(new String(key), new Integer(1));
        } else {
            int v = value.intValue();
            modelCount.put(key, new Integer(++v));
        }
    }

    public void downsize(String key) {
        Integer value = (Integer) modelCount.get(key);
        int v = value.intValue();
        v--;
        if( v == 0 ) {
            modelCount.remove(key);
            savedModels.remove(key);
        } else {
            modelCount.put(key, new Integer(v));
        }
    }

    public HashMap getSavedModels() {
        return savedModels;
    }

    public HashMap getModelCount() {
        return modelCount;
    }

    public HashMap getModelKeys() {
        return modelKeys;
    }

    public String getModekKey(Integer i) {
        return modelKeys.get(i);
    }

    public HeteroCombination getBestModel(Integer i) {
        String key = modelKeys.get(i);
        return savedModels.get(key);
    }

    public int size() {
        return modelKeys.size();
    }

    public void saveModel(HashMap bestmodel) {
        Set KeySet = new HashSet(bestmodel.values());
        HashSet Keys = new HashSet(modelCount.keySet());
        for (Iterator e = KeySet.iterator(); e.hasNext(); ) {
            String key = (String) e.next();
            Keys.remove(key);
        }
        for (Iterator e = Keys.iterator(); e.hasNext(); ) {
            String key = (String) e.next();
            modelCount.remove(key);
            savedModels.remove(key);
        }
        Set KeySet1 = new HashSet(bestmodel.keySet());
        for (Iterator e = KeySet1.iterator(); e.hasNext(); ) {
            Integer key = (Integer) e.next();
            String v = (String) bestmodel.get(key);
            modelKeys.put(key, v);
        }
    }
}
