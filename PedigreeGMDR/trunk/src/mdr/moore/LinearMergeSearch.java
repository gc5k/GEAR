package mdr.moore;

import mdr.Combination;
import mdr.Cell;
import mdr.OneTraitResult;
import mdr.Suite;
import mdr.SavedModels;
import mdr.OneCVSet;
import mdr.BestKFoldCVResult;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.HashMap;
import java.util.List;
import java.util.Set;

import algorithm.CombinationGenerator;
import mdr.data.DataFile;
import mdr.DataFileException;
import algorithm.Subdivision;

import publicAccess.PublicData;
import publicAccess.ToolKit;
import publicAccess.ToolKitException;

/**
 *
 * @author Guo-Bo Chen
 */
public class LinearMergeSearch extends AbstractMergeSearch {

    private int[] SNPIndex;
    private Combination model;
    protected int count;

    public LinearMergeSearch(DataFile dr, Subdivision sd, CombinationGenerator cg, int traits, double[] os, boolean ismooremdr) {
        super(dr, sd, cg, traits, os, ismooremdr);
    }

    public void search(int or) {
        if (testedInteraction.contains(new Integer(or))) {
            return;
        }
        testedInteraction.add(new Integer(or));
        order = or;
        bestKFoldResult = new BestKFoldCVResult(or, numTraits, subdivision.getInterval());
        bestKFoldResultMap.put(new Integer(or), bestKFoldResult);
        count = 0;
        for (int i = 0; i < currBestStatistic.length; i++) {
            for (int j = 0; j < currBestStatistic[i].length; j++) {
                currBestStatistic[i][j] = 0;
            }
        }
        List com = (List) comGenerator.get(or);
        savedModels = new SavedModels();
        for (Iterator e = com.iterator(); e.hasNext();) {
            for (int i = 0; i < cvTestingSet.size(); i++) {
                Combination testingModel = (Combination) cvTestingSet.get(i);
                testingModel.clear();
            }
            model = new Combination();
            String modelName = (String) e.next();
            SNPIndex = ToolKit.StringToIntArray(modelName);
            String c = new String();
            mergeSearch((ArrayList) data.getSample(), c, 0);
            calculate(modelName);
            count++;
        }
        bestSavedModelsMap.put(new Integer(or), savedModels);
        model_summary();
    }

    public void mergeSearch(ArrayList subjects, String combination, int idxMarker) {
        if (idxMarker < SNPIndex.length) {
            HashMap subsets = new HashMap();
            for (Iterator i = subjects.iterator(); i.hasNext();) {
                DataFile.Subject sub = (DataFile.Subject) i.next();
                String m = (String) sub.getGenotype(SNPIndex[idxMarker]);
                if (m.compareTo(PublicData.MissingGenotype) == 0) {
                    continue;
                }
                ArrayList subset = (ArrayList) subsets.get(m);
                if (subset == null) {
                    subset = new ArrayList();
                    subsets.put(m, subset);
                }
                subset.add(sub);
            }
            Set keys = subsets.keySet();
            for (Iterator i = keys.iterator(); i.hasNext();) {
                String key = (String) i.next();
                String com = new String();
                if (combination.length() == 0) {
                    com = key;
                } else {
                    com = combination + PublicData.seperator + key;
                }
                mergeSearch((ArrayList) subsets.get(key), com, idxMarker + 1);
            }
        } // if we've processed all attributes
        else {
            run(combination, subjects);
        }
    }

    protected void run(String com, ArrayList subsample) {
        Suite s = new Suite(subsample, numTraits);
        for (int i = 0; i < numTraits; i++) {
            s.summarize(i, offset[i]);
        }
        model.put(com, s);
        assignKFold(com, subsample);
    }

    protected void assignKFold(String key, ArrayList subsample) {
        for (Iterator e1 = subsample.iterator(); e1.hasNext();) {
            DataFile.Subject sub = (DataFile.Subject) e1.next();
            Integer ID = sub.getIntegerID();
            int d = ((Integer) dataPartitionMap.get(ID)).intValue();
            Combination suiteMap = (Combination) cvTestingSet.get(d);
            Suite S = (Suite) suiteMap.get(key);
            if (S == null) {
                S = new Suite(numTraits);
                suiteMap.put(key, S);
            }
            S.add(sub);
        }

        for (int i = 0; i < numTraits; i++) {
            for (int j = 0; j < subdivision.getInterval(); j++) {
                Combination testingModels = (Combination) cvTestingSet.get(j);
                if (testingModels.containsKey(key)) {
                    Suite testingSuite = (Suite) testingModels.get(key);
                    testingSuite.summarize(i, offset[i]);
                }
            }
        }
    }

    public void calculate(String modelName) {
        Set cellKeys = model.keySet();
        for (int i = 0; i < numTraits; i++) {
            double thresholdRatio = 0;
            try {
                thresholdRatio = data.getTraitRatio(i, SNPIndex);
            } catch (DataFileException E) {
                E.printStackTrace(System.err);
            }
            ArrayList cvSetHolder = new ArrayList();
            boolean flag = false;
            for (int j = 0; j < subdivision.getInterval(); j++) {
                OneCVSet cvSet = new OneCVSet(j, modelName);
                Combination testingModels = (Combination) cvTestingSet.get(j);
                int tr_status;
                int t_status;
                Cell trCell;
                Cell tCell;
                double[] trStatus = new double[cellKeys.size()];
                double[] tStatus = new double[cellKeys.size()];
                int idx = 0;
                for (Iterator e = cellKeys.iterator(); e.hasNext();) {
                    String cellKey = (String) e.next();
                    Suite fullSuite = (Suite) model.get(cellKey);
                    int fullposSubs = fullSuite.getPositiveSubjects(i);
                    int fullnegSubs = fullSuite.getNegativeSubjects(i);
                    double fullposScr = fullSuite.getPositiveScore(i);
                    double fullnegScr = fullSuite.getNegativeScore(i);
                    if (testingModels.containsKey(cellKey)) {
                        Suite testingSuite = (Suite) testingModels.get(cellKey);
                        int pos_Subs = testingSuite.getPositiveSubjects(i);
                        int neg_Subs = testingSuite.getNegativeSubjects(i);
                        double pos_Scr = testingSuite.getPositiveScore(i);
                        double neg_Scr = testingSuite.getNegativeScore(i);
                        tr_status = ToolKit.Acertainment(fullposScr - pos_Scr, fullnegScr - neg_Scr, thresholdRatio);
                        t_status = ToolKit.Acertainment(pos_Scr, neg_Scr, thresholdRatio);
                        trCell = new Cell(fullposSubs - pos_Subs, fullnegSubs - neg_Subs, fullposScr - pos_Scr, fullnegScr - neg_Scr, tr_status);
                        tCell = new Cell(pos_Subs, neg_Subs, pos_Scr, neg_Scr, tr_status);
                    } else {
                        tr_status = ToolKit.Acertainment(fullposScr, fullnegScr, thresholdRatio);
                        t_status = 1 - tr_status;
                        trCell = new Cell(fullposSubs, fullnegSubs, fullposScr, fullnegScr, tr_status);
                        tCell = new Cell(0, 0, 0, 0, -1);
                    }
                    cvSet.addTrainingModel(cellKey, trCell);
                    cvSet.addTestingModel(cellKey, tCell);
                    trStatus[idx] = tr_status;
                    tStatus[idx] = t_status;
                    idx++;
                }
                double cor = ToolKit.getPMCC(trStatus, tStatus);
                double trAccu = 0;
                try {
                    trAccu = ToolKit.Accuracy(cvSet.getTrainingSubdivision());
                } catch (ToolKitException E) {
                    E.printStackTrace(System.err);
                }
                cvSet.setStatistic(PublicData.TrainingAccuIdx, trAccu);
//                cvSet.setStatistic(PublicData.Correlation_Training_Testing, cor);
                if (count == 0) {
                    currBestStatistic[i][j] = trAccu;
                    bestKFoldResult.set(i, j, cvSet);
                    savedModels.save(modelName, model);
                } else if (currBestStatistic[i][j] < trAccu) {
                    currBestStatistic[i][j] = trAccu;
                    String keyAtIJ = bestKFoldResult.getKeyAt(i, j);
                    bestKFoldResult.set(i, j, cvSet);
                    savedModels.save(modelName, model);
                    savedModels.downsize(keyAtIJ);
                }
            }
        }
    }

    public double[][] calculate1(String modelName) {
        Set cellKeys = model.keySet();
        double[][] mean = new double[numTraits][PublicData.NumOfStatistics]; 
        for (int i = 0; i < numTraits; i++) {
            double thresholdRatio = 0;
            try {
                thresholdRatio = data.getTraitRatio(i, SNPIndex);
            } catch (DataFileException E) {
                E.printStackTrace(System.err);
            }
            ArrayList cvSetHolder = new ArrayList();
            boolean flag = false;

            for (int j = 0; j < subdivision.getInterval(); j++) {
                OneCVSet cvSet = new OneCVSet(j, modelName);
                Combination testingModels = (Combination) cvTestingSet.get(j);
                int tr_status;
                int t_status;
                Cell trCell;
                Cell tCell;
                double[] trStatus = new double[cellKeys.size()];
                double[] tStatus = new double[cellKeys.size()];
                int idx = 0;
                for (Iterator e = cellKeys.iterator(); e.hasNext();) {
                    String cellKey = (String) e.next();
                    Suite fullSuite = (Suite) model.get(cellKey);
                    int fullposSubs = fullSuite.getPositiveSubjects(i);
                    int fullnegSubs = fullSuite.getNegativeSubjects(i);
                    double fullposScr = fullSuite.getPositiveScore(i);
                    double fullnegScr = fullSuite.getNegativeScore(i);
                    if (testingModels.containsKey(cellKey)) {
                        Suite testingSuite = (Suite) testingModels.get(cellKey);
                        int pos_Subs = testingSuite.getPositiveSubjects(i);
                        int neg_Subs = testingSuite.getNegativeSubjects(i);
                        double pos_Scr = testingSuite.getPositiveScore(i);
                        double neg_Scr = testingSuite.getNegativeScore(i);
                        tr_status = ToolKit.Acertainment(fullposScr - pos_Scr, fullnegScr - neg_Scr, thresholdRatio);
                        t_status = ToolKit.Acertainment(pos_Scr, neg_Scr, thresholdRatio);
                        trCell = new Cell(fullposSubs - pos_Subs, fullnegSubs - neg_Subs, fullposScr - pos_Scr, fullnegScr - neg_Scr, tr_status);
                        tCell = new Cell(pos_Subs, neg_Subs, pos_Scr, neg_Scr, tr_status);
                    } else {
                        tr_status = ToolKit.Acertainment(fullposScr, fullnegScr, thresholdRatio);
                        t_status = 1 - tr_status;
                        trCell = new Cell(fullposSubs, fullnegSubs, fullposScr, fullnegScr, tr_status);
                        tCell = new Cell(0, 0, 0, 0, -1);
                    }
                    cvSet.addTrainingModel(cellKey, trCell);
                    cvSet.addTestingModel(cellKey, tCell);
                    trStatus[idx] = tr_status;
                    tStatus[idx] = t_status;
                    idx++;
                }
                double cor = ToolKit.getPMCC(trStatus, tStatus);
                double trAccu = 0;
                double tAccu = 0;
                try {
                    trAccu = ToolKit.Accuracy(cvSet.getTrainingSubdivision());
                    mean[i][PublicData.TrainingAccuIdx] += trAccu;
                    tAccu = ToolKit.Accuracy(cvSet.getTestingSubdivision());
                    mean[i][PublicData.TestingAccuIdx] += tAccu;
                } catch (ToolKitException E) {
                    E.printStackTrace(System.err);
                }
                cvSet.setStatistic(PublicData.TrainingAccuIdx, trAccu);
                cvSet.setStatistic(PublicData.TestingAccuIdx, tAccu);
            }
            mean[i][PublicData.TrainingAccuIdx] /= subdivision.getInterval();
            mean[i][PublicData.TestingAccuIdx] /= subdivision.getInterval();
        }
        return mean;
    }

    private void model_summary() {
        HashMap bestModelKeyMap = new HashMap();
        for (int i = 0; i < numTraits; i++) {
            for (int j = 0; j < subdivision.getInterval(); j++) {
                OneCVSet cvSet = bestKFoldResult.get(i, j);
                double testingAccu = 0;
                try {
                    testingAccu = ToolKit.Accuracy(cvSet.getTestingSubdivision());
                } catch (ToolKitException E) {
                    E.printStackTrace(System.err);
                }
                cvSet.setStatistic(PublicData.TestingAccuIdx, testingAccu);
//                System.out.println(cvSet.getCombination() + " " + cvSet.getStatistic(PublicData.TestingAccuIdx) + " " + cvSet.getStatistic(PublicData.TrainingAccuIdx));
            }
            OneTraitResult i_thResult = (OneTraitResult) bestKFoldResult.get(i);
            i_thResult.summarise();
            bestModelKeyMap.put(new Integer(i), i_thResult.getBestModelKey());
            stats[i] = i_thResult.getStatistic();
        }
        SavedModels sModels = (SavedModels) bestSavedModelsMap.get(new Integer(order));
        sModels.saveModel(bestModelKeyMap);
    }

    public double[][] singleBest(String com) {
        for (int i = 0; i < cvTestingSet.size(); i++) {
            Combination testingModel = (Combination) cvTestingSet.get(i);
            testingModel.clear();
        }
        model = new Combination();
        SNPIndex = ToolKit.StringToIntArray(com);
        String c = new String();
        mergeSearch((ArrayList) data.getSample(), c, 0);
        return calculate1(com);
    }

    public void summarise() {
        HashMap bestModelKeyMap = new HashMap();
        for (int i = 0; i < numTraits; i++) {
            for (int j = 0; j < subdivision.getInterval(); j++) {
                OneCVSet cvSet = bestKFoldResult.get(i, j);
            }
            OneTraitResult i_thResult = (OneTraitResult) bestKFoldResult.get(i);
            System.out.println("Trait name: " + data.getTraitName(i));
            System.out.println(i_thResult);
        }
    }
}
