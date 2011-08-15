package family.mdr.result;

import java.util.HashMap;

import family.mdr.MDRConstant;
import util.NewIt;

/**
 *
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class OneCVSet {

    int cvIdx;
    String model;
    double[] statistic;
    double cor;
    HashMap<String, Cell> trainingMap;
    HashMap<String, Cell> testingMap;

    public OneCVSet(int i, String com) {
        cvIdx = i;
        model = com;
        statistic = new double[MDRConstant.NumStats];
        trainingMap = NewIt.newHashMap();
        testingMap = NewIt.newHashMap();
    }

    public OneCVSet(HashMap<String, Cell> trn, HashMap<String, Cell> tst) {
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

    public HashMap<String, Cell> getTrainingSubdivision() {
        return trainingMap;
    }

    public HashMap<String, Cell> getTestingSubdivision() {
        return testingMap;
    }

    public int getCVIndex() {
        return cvIdx;
    }

    public String getModel() {
        return model;
    }
}
