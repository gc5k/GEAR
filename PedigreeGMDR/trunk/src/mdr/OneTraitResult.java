package mdr;

import mdr.moore.*;
import mdr.OneCVSet;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import mdr.Cell;

import publicAccess.PublicData;
/**
 *
 * @author Guo-Bo Chen
 */
public class OneTraitResult extends ArrayList {

    private int traitIdx;
    private double[] statistic;
    private double[] bestModelStatistic;
    private HashMap modelCount = new HashMap();
    private String bestModelKey;
    private int cvConsistency;
    private StringBuffer sb;
    public OneTraitResult(int i) {
        traitIdx = i;
    }

    public void print() {
        int cv = 0;
        for (Iterator e = iterator(); e.hasNext();) {
            OneCVSet cvPair = (OneCVSet) e.next();
            System.out.println("CV result " + cv + ", Testing Accuracy " + cvPair.getStatistic(PublicData.TestingAccuIdx) 
                    + ", Training Accuracy " + cvPair.getStatistic(PublicData.TrainingAccuIdx));
            HashMap trMap = cvPair.getTrainingSubdivision();
            TreeSet trkeyset = new TreeSet(trMap.keySet());
            for (Iterator se = trkeyset.iterator(); se.hasNext();) {
                String key = (String) se.next();
                    System.out.println(key);
                Cell trCell = (Cell) trMap.get(key);
                    System.out.println(trCell);
            }
        }
        cv++;
    }

    public int getTraitIndex() {
        return traitIdx;
    }

    public HashMap getModelCount() {
        return modelCount;
    }

    public double[] getStatistic() {
        return statistic;
    }

    public double[] getBestModelStatistic() {
    	return bestModelStatistic;
    }

    public String getBestModelKey() {
        return bestModelKey;
    }

    public void summarise() {
        statistic = new double[PublicData.NumOfStatistics];
        bestModelStatistic = new double[PublicData.NumOfStatistics];

        for (int i = 0; i < size(); i++) {
            OneCVSet cvSet = (OneCVSet) get(i);
            for(int j = 0; j < PublicData.NumOfStatistics; j++ ) {
                statistic[j] += cvSet.getStatistic(j);
            }

            String key = cvSet.getCombination();
            Integer count = (Integer) modelCount.get(key);
            if (count == null) {
                modelCount.put(new String(key), new Integer(1));
            } else {
                int v = count.intValue();
                modelCount.put(key, new Integer(++v));
            }
        }

        for (int i = 0; i < PublicData.NumOfStatistics; i++) {
            statistic[i] /= size();
        }
        
        int big = 0;
        Set keys = modelCount.keySet();
        for (Iterator e = keys.iterator(); e.hasNext(); ) {
            Integer count = (Integer) modelCount.get((String) e.next());
            if(big < count.intValue()) {
                big = count.intValue();
            }
        }
        HashMap majorModel = new HashMap();
        ArrayList majorKey = new ArrayList();
        for (Iterator e = keys.iterator(); e.hasNext(); ) {
            String key = (String) e.next();
            Integer count = (Integer) modelCount.get(key);
            if(big == count.intValue()) {
                majorModel.put(key, count);
                majorKey.add(key);
            }
        }

        double[][] _bestModelStats = new double[majorKey.size()][PublicData.NumOfStatistics];
        for (int i = 0; i < size(); i++) {
            OneCVSet cvSet = (OneCVSet) get(i);
            String key = cvSet.getCombination();
            if (!majorKey.contains(key)) {
                continue;
            }
            int idx = majorKey.indexOf(key);
            for (int j = 0; j < PublicData.NumOfStatistics; j++) {
            	_bestModelStats[idx][j] += cvSet.getStatistic(j);
            }
        }
        double[] bigStats = new double[PublicData.NumOfStatistics];
        int idx = 0;
        for (Iterator e = majorKey.iterator(); e.hasNext();) {
            String key = (String) e.next();
            for (int j = 0; j < _bestModelStats[idx].length; j++) {
            	_bestModelStats[idx][j] /= ((Integer) (majorModel.get(key))).intValue();
                if (bigStats[j] < _bestModelStats[idx][j]) {
                    bigStats[j] = _bestModelStats[idx][j];
                }
            }
            idx++;
        }
        for (int i = 0; i < PublicData.NumOfStatistics; i++) {
            int c = 0;
            for (int j = 0; j < _bestModelStats.length; j++) {
                if ((bigStats[i] - _bestModelStats[j][i]) < PublicData.epsilon) {
                    c++;
                    bestModelKey = (String) majorKey.get(j);
                }
            }
            if (c == 1) {
                break;
            }
        }
        int bestidx = majorKey.indexOf(bestModelKey);
        cvConsistency = big;
        System.arraycopy(_bestModelStats[bestidx], 0, bestModelStatistic, 0, _bestModelStats[bestidx].length);
    }

    public String toString() {
    	sb = new StringBuffer();
    	sb.append("Statistics of the CV model for trait " + traitIdx + System.getProperty("line.separator"));
        for (int i = 0; i < size(); i++) {
            OneCVSet cvSet = (OneCVSet) get(i);
            for(int j = 0; j < PublicData.NumOfStatistics; j++ ) {
                statistic[j] += cvSet.getStatistic(j);
            }

            String key = cvSet.getCombination();
            Integer count = (Integer) modelCount.get(key);
            if (count == null) {
                modelCount.put(new String(key), new Integer(1));
            } else {
                int v = count.intValue();
                modelCount.put(key, new Integer(++v));
            }
        }
        sb.append("Statistics of the best model" + System.getProperty("line.separator"));
        sb.append("Best model: " + bestModelKey + System.getProperty("line.separator"));
        sb.append("Testing Accuracy: " + bestModelStatistic[PublicData.TestingAccuIdx] + System.getProperty("line.separator"));
        sb.append("Training Accuracy: " + bestModelStatistic[PublicData.TrainingAccuIdx] + System.getProperty("line.separator"));
        sb.append("Cross-validation consistency: " + cvConsistency + "/" + size());

        return sb.toString();
    }
}
