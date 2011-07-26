
package mdr.algorithm;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Random;

import util.NewIt;

/**
 *
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class Subdivision {
    private HashMap<Integer, Integer> divisionMap;    // a map from subject IDs to interval IDs.
    private int populationSize; //0 for DataFile, 1 for IMPopulation
    private long seed;
    private int interval;                           // number of intervals

    public Subdivision(int intl, long sd, int ps) {
        seed = sd;
        interval = intl;
        populationSize = ps;
        divisionMap = NewIt.newHashMap();
    }

    public void RandomPartition() {
        ArrayList<Integer> ID = NewIt.newArrayList();
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

    public HashMap<Integer, Integer> getDivision() {
        return divisionMap;
    }
    
    public int getInterval() {
        return interval;
    }
}
