package mdrAlternative;

import mdr.Combination;
import mdr.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Iterator;
import java.util.HashMap;
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
public class DynamicSearch extends AbstractSearch {

    public DynamicSearch(DataFile dr, CombinationGenerator cg, int[] sI) {
        super(dr, cg, sI);
    }

    public void search(int or) {
        if( modelMap.size() == 0) {
            for (int i = 1; i <= or; i++) {
                stretch(i);
            }
        } else if(modelMap.containsKey(new Integer(or))) {
            return;
        } else {
            int begin = 1;
            for ( int i = or; i>=1; i--) {
                if ( modelMap.containsKey(new Integer(or))) {
                    begin = i;
                    break;
                }
            }
            for (int i = begin + 1; i<= or; i++) {
                stretch(i);
            }
        }
    }

    private void stretch(int or) {
        List subjects = (ArrayList) data.getSample();
        Map localModel = new HashMap();
        if (or == 1) {
            for (int i = 0; i < data.getMarkerNum(); i++) {
                Map subsets = new HashMap();
                for (Iterator e = subjects.iterator(); e.hasNext();) {
                    DataFile.Subject sub = (DataFile.Subject) e.next();
                    String m = (String) sub.getGenotype(i);
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
                Combination model = new Combination();
                Set keys = subsets.keySet();
                for (Iterator e = keys.iterator(); e.hasNext();) {
                    String key = (String) e.next();
                    ArrayList subsample = (ArrayList) subsets.get(key);
                    Suite s = new Suite(subsample, scrIdx.length);
                    model.put(key, s);
                }
                localModel.put((new Integer(i)).toString(), model);
            }
        } else {
            Integer mkey = new Integer(or - 1);
            Map models = (HashMap) modelMap.get(mkey);
            Set keys = models.keySet();
            for (Iterator e = keys.iterator(); e.hasNext();) {
                String key = (String) e.next();
                int[] idx = ToolKit.StringToIntArray(key);
                for (int i = idx[idx.length - 1] + 1; i < data.getMarkerNum(); i++) {
                    Combination tmpModel = (Combination) models.get(key);
                    Set genosets = tmpModel.keySet();
                    Combination model = new Combination();
                    for (Iterator e1 = genosets.iterator(); e1.hasNext();) {
                        String geno = (String) e1.next();
                        Suite s = (Suite) tmpModel.get(geno);
                        Map submodel = StretchSuite(geno, i, s);
                        Set subgenosets = submodel.keySet();
                        for (Iterator e2 = subgenosets.iterator(); e2.hasNext();) {
                            String subgeno = (String) e2.next();
                            Suite S = (Suite) submodel.get(subgeno);
                            model.put(subgeno, S);
                        }
                    }
                    String newkey = key + PublicData.seperator + Integer.toString(i);
                    localModel.put(newkey, model);
                }
            }
        }
        modelMap.put(new Integer(or), localModel);
    }

    public void existOrder(int or) {
        Map hm = null;
        if (modelMap.size() >= or) {
            hm = (HashMap) modelMap.get(new Integer(or));
        }
    }

    public Map StretchSuite(String key, int idxMarker, Suite s) {
        Map subsets = new HashMap();
        List subs = s.getSubjects();
        for (Iterator e = subs.iterator(); e.hasNext();) {
            DataFile.Subject sub = (DataFile.Subject) e.next();
            String currGeno = (String) sub.getGenotype(idxMarker);
            if(currGeno.compareTo(PublicData.MissingGenotype) == 0) {
                continue;
            }
            String geno = key + PublicData.seperator + currGeno;
            List subset = (ArrayList) subsets.get(geno);
            if (subset == null) {
                subset = new ArrayList();
                subsets.put(geno, subset);
            }
            subset.add(sub);
        }
        Map suiteMap = new HashMap();
        Set genos = subsets.keySet();
        for (Iterator e = genos.iterator(); e.hasNext();) {
            String geno = (String) e.next();
            ArrayList subsample = (ArrayList) subsets.get(geno);
            Suite S = new Suite(subsample, scrIdx.length);
            suiteMap.put(geno, S);
        }
        return suiteMap;
    }
}
