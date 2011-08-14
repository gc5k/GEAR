package family.mdr.result;

import java.util.HashMap;

import family.mdr.result.Combination;


import util.NewIt;

/**
 *
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class SavedModels {

    HashMap<String, Combination> savedModels = NewIt.newHashMap();
    HashMap<String, Integer> modelCount = NewIt.newHashMap();
    HashMap<Integer, String> modelKeys = NewIt.newHashMap();//keep the key of the best model to each trait.

    public SavedModels() {}

    public void save(String key, Combination model, String deleteModel) {
        if( !savedModels.containsKey(key)) {
            savedModels.put(key, model);
        }
        Integer value = modelCount.get(key);
        if( value == null ) {
            modelCount.put(key, new Integer(1));
        } else {
            value++;
            modelCount.put(key, value);
        }
        if(deleteModel != null) downsize(deleteModel);
    }

    private void downsize(String key) {
        Integer value = (Integer) modelCount.get(key);
        value--;;
        if( value == 0 ) {
            modelCount.remove(key);
            savedModels.remove(key);
        } else {
            modelCount.put(key, value);
        }
    }
    
    public HashMap<String, Combination> getSavedModels() {
        return savedModels;
    }
    
    public HashMap<String, Integer> getModelCount() {
        return modelCount;
    }
    
    public HashMap<Integer, String> getModelKeys() {
        return modelKeys;
    }

    public String getModelKey(Integer i) {
        return modelKeys.get(i);
    }

    public Combination getBestModel(Integer i) {
        String key = modelKeys.get(i);
        return savedModels.get(key);
    }

    public int size() {
        return modelKeys.size();
    }

    public void saveModel(HashMap<Integer, String> bestmodel) {

        for (String key : bestmodel.values()) {
            modelCount.remove(key);
        }

        for (String key : modelCount.keySet() ) {
            modelCount.remove(key);
            savedModels.remove(key);
        }

        for (Integer key : bestmodel.keySet()) {
            String v = bestmodel.get(key);
            modelKeys.put(key, v);
        }
    }
}
