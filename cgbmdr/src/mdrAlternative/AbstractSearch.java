package mdrAlternative;

import java.util.AbstractMap;
import java.util.Collections;
import java.util.Iterator;
import java.util.Map;
import java.util.HashMap;
import java.util.Set;
import java.util.TreeSet;

import algorithm.CombinationGenerator;
import mdr.data.DataFile;
import algorithm.Subdivision;
import mdr.Suite;

/**
 *
 * @author Guo-Bo Chen
 */
public abstract class AbstractSearch extends AbstractMap {

    protected HashMap modelMap = new HashMap();

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

    public Set entrySet() {
        return Collections.unmodifiableSet(modelMap.entrySet());
    }

    public HashMap getModelMap(Integer or) {
        return (HashMap) modelMap.get(or);
    }

    public HashMap getModelMap() {
        return modelMap;
    }
    
    public int size() {
        return modelMap.size();
    }
    
    public void testPrint() {
        Set keys = modelMap.keySet();
        for (Iterator e = keys.iterator(); e.hasNext();) {
            Integer key = (Integer) e.next();
            System.out.println("Order of interaction " + key);
            Map hm = (HashMap) modelMap.get(key);
            Set modelkeys = new TreeSet(hm.keySet());
            for (Iterator e1 = modelkeys.iterator(); e1.hasNext();) {
                String modelkey = (String) e1.next();
                System.out.println("model " + modelkey);
                Map model = (HashMap) hm.get(modelkey);
                Set genoset = new TreeSet(model.keySet());
                for (Iterator e2 = genoset.iterator(); e2.hasNext();) {
                    String geno = (String) e2.next();
                    Suite s = (Suite) model.get(geno);
                    System.out.println(geno);
                    System.out.println(s);
                }
            }
        }
    }
}
