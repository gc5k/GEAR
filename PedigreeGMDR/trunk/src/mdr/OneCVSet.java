package mdr;

import mdr.moore.*;
import java.util.HashMap;

import mdr.Cell;

import publicAccess.PublicData;
/**
 *
 * @author Guo-Bo Chen
 */
public class OneCVSet {

    int cvIdx;
    String combination;
    double[] statistic;
    double cor;
    HashMap trainingMap;
    HashMap testingMap;

    public OneCVSet(int i, String com) {
        cvIdx = i;
        combination = com;
        statistic = new double[PublicData.NumOfStatistics];
        trainingMap = new HashMap();
        testingMap = new HashMap();
    }

    public OneCVSet(HashMap trn, HashMap tst) {
        trainingMap = trn;
        testingMap = tst;
    }

    public double getStatistic(int idx) {
        return statistic[idx];
    }
    
    public void setStatistic(int idx, double value) {
        statistic[idx] = value;
    }

    public void addTrainingModel(String key, Cell c) {
        trainingMap.put(key, c);
    }

    public void addTestingModel(String key, Cell c) {
        testingMap.put(key, c);
    }

    public HashMap getTrainingSubdivision() {
        return trainingMap;
    }

    public HashMap getTestingSubdivision() {
        return testingMap;
    }

    public int getCVIndex() {
        return cvIdx;
    }

    public String getCombination() {
        return combination;
    }
}
