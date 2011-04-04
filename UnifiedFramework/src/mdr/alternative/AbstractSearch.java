package mdr.alternative;

import java.util.AbstractMap;
import java.util.Collections;
import java.util.Map;
import java.util.HashMap;
import java.util.Set;
import java.util.TreeSet;

import algorithm.CombinationGenerator;
import mdr.data.DataFile;
import mdr.Suite;
import mdr.Combination;
import util.NewIt;
/**
 *
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public abstract class AbstractSearch extends AbstractMap<Integer, HashMap<String, Combination>> {

    protected HashMap<Integer, HashMap<String, Combination>> modelMap = NewIt.newHashMap();

    protected int order;
    protected int[] scrIdx;
    protected CombinationGenerator comGenerator;
    protected DataFile data;

    public AbstractSearch(DataFile dr, CombinationGenerator cg, int[] sI) {
        data = dr;
        comGenerator = cg;
        scrIdx = new int[sI.length];
        System.arraycopy(sI, 0, scrIdx, 0, sI.length);
    }

    public abstract void search(int or);

    public Set<Map.Entry<Integer, HashMap<String, Combination>>> entrySet() {
        return Collections.unmodifiableSet( modelMap.entrySet());
    }

    public HashMap<String, Combination> getModelMap(Integer or) {
        return (HashMap<String, Combination>) modelMap.get(or);
    }

    public HashMap<Integer, HashMap<String, Combination>> getModelMap() {
        return modelMap;
    }
    
    public int size() {
        return modelMap.size();
    }
    
    public void testPrint() {
//        for (Iterator e = keys.iterator(); e.hasNext();) {
        for(Integer key:modelMap.keySet()) {
            System.out.println("Order of interaction " + key);
            Map<String, Combination> hm = modelMap.get(key);
            Set<String> modelkeys = new TreeSet<String>(hm.keySet());
//            for (Iterator e1 = modelkeys.iterator(); e1.hasNext();) {
            for(String modelkey:modelkeys) {
                System.out.println("model " + modelkey);
                Map<String, Suite> model = hm.get(modelkey);
                Set<String> genoset = new TreeSet<String>(model.keySet());
                for (String geno:genoset) {
                    Suite s = (Suite) model.get(geno);
                    System.out.println(geno);
                    System.out.println(s);
                }
            }
        }
    }
}
