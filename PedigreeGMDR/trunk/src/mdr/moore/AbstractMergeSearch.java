package mdr.moore;

import mdr.Combination;
import mdr.SavedModels;
import mdr.BestKFoldCVResult;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import publicAccess.ToolKit;

import algorithm.CombinationGenerator;
import mdr.data.DataFile;
import algorithm.Subdivision;

/**
 *
 * @author Guo-Bo Chen
 */
public abstract class AbstractMergeSearch {

    protected boolean isMooreMDR;
    protected int order;
    protected int numTraits;
    protected double[] offset;
    protected double[][] currBestStatistic;//keeps currently the best statistic given the order of interaction

    protected HashSet<Integer> testedInteraction = new HashSet();//Interaction has already tested
    protected HashMap dataPartitionMap;
    protected HashMap<Integer, BestKFoldCVResult> bestKFoldResultMap = new HashMap();//keeps the best results under each interaction
    protected HashMap<Integer, SavedModels> bestSavedModelsMap = new HashMap();//keeps the best model under each interaction
    protected ArrayList<Combination> cvTestingSet = new ArrayList();//always keeps the K testing models of the current combination

    protected BestKFoldCVResult bestKFoldResult;
    protected CombinationGenerator comGenerator;
    protected DataFile data;
    protected SavedModels savedModels;
    protected Subdivision subdivision;

    public AbstractMergeSearch(DataFile dr, Subdivision sd, CombinationGenerator cg, int traits, double[] os, boolean ismooremdr) {
        data = dr;
        subdivision = sd;
        comGenerator = cg;
        numTraits = traits;
        offset = new double[os.length];
        System.arraycopy(os, 0, offset, 0, os.length);
        isMooreMDR = ismooremdr;

        for (int i = 0; i < subdivision.getInterval(); i++) {
            Combination testingMap = new Combination();
            cvTestingSet.add(testingMap);
        }

        dataPartitionMap = subdivision.getDivision();
        currBestStatistic = new double[numTraits][];
        for (int i = 0; i < numTraits; i++) {
            currBestStatistic[i] = new double[subdivision.getInterval()];
        }
    }

    public int size() {
        return testedInteraction.size();
    }

    public boolean isMooreMDR() {
        return isMooreMDR;
    }
    
    public SavedModels getBestSavedModelAtOrder(Integer i) {
        return bestSavedModelsMap.get(i);
    }
    
    public int[] getBestModel(int o, int trIdx) {
    	SavedModels sm = bestSavedModelsMap.get(new Integer(o));
    	String m = sm.getModekKey(trIdx);
    	return ToolKit.StringToIntArray(m);
    }

    public abstract void search(int or);
    protected abstract void assignKFold(String key, ArrayList subsample);
    public abstract void summarise();
}
