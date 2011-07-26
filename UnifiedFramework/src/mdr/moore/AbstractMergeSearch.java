package mdr.moore;

import java.util.ArrayList;
import java.util.HashMap;


import mdr.algorithm.Subdivision;
import mdr.arsenal.ToolKit;
import mdr.data.DataFile;
import mdr.moore.statistic.MDRStatistic;
import mdr.result.BestKFoldCVResult;
import mdr.result.Combination;
import mdr.result.SavedModels;


import util.NewIt;

/**
 *
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public abstract class AbstractMergeSearch {
	protected int order;
    protected double[] currBestStats;//keeps currently the best statistic

    protected ArrayList<Combination> cvTestingSet = NewIt.newArrayList();//always keeps the K testing models of the current combination
    protected BestKFoldCVResult bestKFold;

    protected HashMap<Integer, Integer> dataPartitionMap;
    protected DataFile data;
    protected SavedModels savedModels;
    protected Subdivision subdivision;

    protected double[] stats;

	protected HashMap<String, MDRStatistic> heteroresult = null;
	
    public AbstractMergeSearch(DataFile dr, Subdivision sd) {
        data = dr;
        subdivision = sd;
        for (int i = 0; i < subdivision.getInterval(); i++) {
            Combination testingMap = new Combination();
            cvTestingSet.add(testingMap);
        }

        dataPartitionMap = subdivision.getDivision();
        currBestStats = new double[subdivision.getInterval()];
    }

    public double[] getStats() {
    	return stats;
    }

    public int[] getBestModel() {
    	String m = bestKFold.getBestModel();
    	return ToolKit.StringToIntArray(m);
    }

    public String getBestModelKey() {
    	return bestKFold.getBestModel();
    }

    public HashMap<String, MDRStatistic> getMDRResult() {
    	return heteroresult;
    }

    public abstract void search(int or, ArrayList<String> modelspace);
    protected abstract void assignKFold(String key, ArrayList<DataFile.Subject> subsample);
    public abstract void summarise();
    public abstract String toString();

}
