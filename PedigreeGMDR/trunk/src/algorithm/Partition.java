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

    ArrayList<HashSet<Integer>> testingList = new ArrayList<HashSet<Integer>>();    // list of partitions
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
            testingList.add(new HashSet<Integer>());
        }
    }

    public void partition(int idx, int method) {
        if (method == PublicData.RandomPartition) {
            randomPartition();
        } else if(method == PublicData.UnpairedPartition ) {
            unpairedPartition(idx);
        }
    }

    private void unpairedPartition(int idx) {
        ArrayList<Integer> affecteds = new ArrayList<Integer>();
        ArrayList<Integer> unaffecteds = new ArrayList<Integer>();
        for (DataFile.Subject sub:data.getSample()) {
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
            ((HashSet<Integer>) testingList.get(i % interval)).add(affecteds.get(i));
        }

        for (int i = 0; i < unaffecteds.size(); ++i) {
            ((HashSet<Integer>) testingList.get((i + affecteds.size()) % interval)).add(unaffecteds.get(i));
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
        testingList.add(idx, (HashSet<Integer>) o);
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
