
package algorithm;

import mdr.data.DataFile;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Random;

import im.population.IMPopulation;
/**
 *
 * @author Guo-Bo Chen
 */
public class Subdivision {
    private HashMap divisionMap;    // a map from subject IDs to interval IDs.
    private DataFile data;
    private IMPopulation imp;
    private int populationSize; //0 for DataFile, 1 for IMPopulation
    private long seed;
    private int interval;                           // number of intervals

    public Subdivision(int intl, long sd, DataFile dr) {
        data = dr;
        seed = sd;
        interval = intl;
        populationSize = data.size();
        divisionMap = new HashMap();
    }

    public Subdivision(int intl, long sd, IMPopulation ip) {
        imp = ip;
        seed = sd;
        interval = intl;
        populationSize = imp.IndividualNumber();
        divisionMap = new HashMap();
    }

    private void UnPairedPartition(int idx, double[] offset) {
        ArrayList affecteds = new ArrayList();
        ArrayList unaffecteds = new ArrayList();
        ArrayList temp_data = data.getSample();
        for (Iterator i = temp_data.iterator(); i.hasNext();) {
            DataFile.Subject sub = (DataFile.Subject) i.next();
            Double subscore = sub.getDoubleScore(idx);
            if (subscore != null) {
                double scr = subscore - offset[idx];
                if ((scr) >= 0) {
                    affecteds.add(sub.getIntegerID());
                } else {
                    unaffecteds.add(sub.getIntegerID());
                }
            }
        }
        Random rnd = new Random(seed);
        if (rnd != null) {
            Collections.shuffle(affecteds, rnd);
            Collections.shuffle(unaffecteds, rnd);
        }

        for (int i = 0; i < affecteds.size(); ++i) {
            Integer value = new Integer(i%interval);
            divisionMap.put((Integer) affecteds.get(i), value);
        }

        for (int i = 0; i < unaffecteds.size(); ++i) {
            Integer value = new Integer((i + affecteds.size())%interval);
            divisionMap.put((Integer) unaffecteds.get(i), value);
        }
    }
    
    public void RandomPartition() {
        ArrayList ID = new ArrayList();
        for (int i = 0; i < populationSize; i++) {
            ID.add(new Integer(i));
        }

        Random rnd = new Random(seed);
        Collections.shuffle(ID, rnd);

        for (int i = 0; i < populationSize; i++) {
            Integer value = new Integer(i%interval);
            divisionMap.put(ID.get(i), value);
        }
    }

    public HashMap getDivision() {
        return divisionMap;
    }
    
    public int getInterval() {
        return interval;
    }
}
