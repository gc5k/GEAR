package mdr.algorithm;

import mdr.MDRConstant;
import mdr.data.DataFile;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Random;

import util.NewIt;

/**
 *
 * @author Guo-Bo Chen
 */
public class Partition {

    ArrayList<HashSet<Integer>> testingList = NewIt.newArrayList();    // list of partitions
    int populationSize;
    int interval;
    long seed;
    int numTraits;
    DataFile data;
    public Partition(int interval, long seed, DataFile ds, int traits) {
        this.interval = interval;
        this.seed = seed;
        populationSize = ds.size();
        data = ds;
          
        numTraits = traits;

        for (int i = 0; i < interval; i++) {
            testingList.add(new HashSet<Integer>());
        }
    }

    public void partition(int idx, int method) {
        if (method == MDRConstant.RandomPartition) {
            randomPartition();
        } else if(method == MDRConstant.UnpairedPartition ) {
            unpairedPartition(idx);
        }
    }

    private void unpairedPartition(int idx) {
        ArrayList<Integer> affecteds = new ArrayList<Integer>();
        ArrayList<Integer> unaffecteds = new ArrayList<Integer>();
        for (DataFile.Subject sub:data.getSample()) {
            double subscore = sub.getSelectedScore();
            if (subscore != Double.NaN) {
                double scr = subscore;
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
        ArrayList<Integer> ID = NewIt.newArrayList();
        for (int i = 0; i < populationSize; i++) {
            ID.add(new Integer(i));
        }

        Random rnd = new Random(seed);
        Collections.shuffle(ID, rnd);

        for (int i = 0; i < populationSize; i++) {
            HashSet<Integer> hs = (HashSet<Integer>) testingList.get(i % interval);
            hs.add((Integer) ID.get(i));
        }
    }

    public void add(int idx, HashSet<Integer> o) {
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
        for (HashSet<Integer> hs : testingList) {
            for (Integer e : hs) {
                System.out.print(e + " ");
            }
            System.out.println();
        }
    }
}
