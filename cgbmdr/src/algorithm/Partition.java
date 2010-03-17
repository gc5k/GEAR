package algorithm;

import mdr.data.DataFile;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Random;

import publicAccess.PublicData;
/**
 *
 * @author Guo-Bo Chen
 */
public class Partition extends AbstractList {

    ArrayList testingList = new ArrayList();    // list of partitions
    int populationSize;
    int interval;
    long seed;
    int numTraits;
    double[] offset;
    DataFile data;
    public Partition(int interval, long seed, DataFile ds, int traits, double[] os) {
        this.interval = interval;
        this.seed = seed;
        populationSize = ds.size();
        data = ds;
          
        numTraits = traits;
        offset = new double[os.length];
        System.arraycopy(os, 0, offset, 0, os.length);

        for (int i = 0; i < interval; i++) {
            testingList.add(new HashSet());
        }
    }

    public void partition(int idx, int method) {
        if (method == PublicData.RandomPartition) {
            randomPartition();
        } else if(method == PublicData.UnpairedPartition ) {
            unpairedPartition(idx);
        } else {
            pairedPartition();
        }
    }

    private void unpairedPartition(int idx) {
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
            ((HashSet) testingList.get(i % interval)).add(affecteds.get(i));
        }

        for (int i = 0; i < unaffecteds.size(); ++i) {
            ((HashSet) testingList.get((i + affecteds.size()) % interval)).add(unaffecteds.get(i));
        }
    }

    private void pairedPartition() {
        ArrayList pairs = new ArrayList((populationSize + 1) / 2);

        for (int i = 0; i < populationSize; i += 2) {
            Object[] pair;
            int j = i + 1;
            if (j < populationSize) {
                pair = new Object[2];
                pair[0] = data.get(i);
                pair[1] = data.get(i + 1);
            } else {
                pair = new Object[1];
                pair[0] = data.get(i);
            }
            
            pairs.add(pair);
        }
        Random rnd = new Random(seed);
        if (rnd != null) {
            Collections.shuffle(pairs, rnd);
        }

        for (int i = 0; i < pairs.size(); ++i) {
            HashSet set = (HashSet) testingList.get(i % interval);
            Object[] pair = (Object[]) pairs.get(i);
            set.add(pair[0]);
            if (pair.length > 1)
                set.add(pair[1]);
        }
    }

    private void randomPartition() {
        ArrayList ID = new ArrayList();
        for (int i = 0; i < populationSize; i++) {
            ID.add(new Integer(i));
        }

        Random rnd = new Random(seed);
        Collections.shuffle(ID, rnd);

        for (int i = 0; i < populationSize; i++) {
            HashSet hs = (HashSet) testingList.get(i % interval);
            hs.add((Integer) ID.get(i));
        }
    }

    public void add(int idx, Object o) {
        modCount++;
        testingList.add(idx, o);
    }

    public Object get(int idx) {
        return testingList.get(idx);
    }

    public Object remove(int idx) {
        return testingList.remove(idx);
    }

    public int size() {
        return testingList.size();
    }

    public void print() {
        for (Iterator e = testingList.iterator(); e.hasNext();) {
            HashSet hs = (HashSet) e.next();
            for (Iterator ee = hs.iterator(); ee.hasNext();) {
                System.out.print((Integer) ee.next() + " ");
            }
            System.out.println();
        }
    }
}
