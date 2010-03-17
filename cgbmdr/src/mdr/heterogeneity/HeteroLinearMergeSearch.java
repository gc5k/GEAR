package mdr.heterogeneity;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.HashMap;
import java.util.List;
import java.util.Set;

import algorithm.CombinationGenerator;
import mdr.data.DataFile;
import mdr.DataFileException;
import algorithm.Subdivision;

import mdr.Cell;
import mdr.SavedModels;
import mdr.Combination;
import mdr.OneCVSet;
import mdr.OneTraitResult;
import mdr.Suite;
import mdr.BestKFoldCVResult;

import publicAccess.PublicData;
import publicAccess.ToolKit;
import publicAccess.ToolKitException;

/**
 *
 * @author Guo-Bo Chen
 */
public class HeteroLinearMergeSearch extends AbstractHeteroMergeSearch {

    protected DataFile data;

    public HeteroLinearMergeSearch(DataFile dr, Subdivision sd, CombinationGenerator cg, int traits, double[] os, boolean ismooremdr) {
        super(sd, cg, traits, os, ismooremdr);
        data = dr;
    }

    /**
     *
     * @param or    number of loci
     */
    public void search(int or, int pheIdx) {
        if (testedInteraction.contains(new Integer(or))) {
            return;
        }
        testedInteraction.add(new Integer(or));
        bestKFoldResult = new BestKFoldCVResult(or, numTraits, subdivision.getInterval());
        bestKFoldResultMap.put(new Integer(or), bestKFoldResult);
        order = or;
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
            mergeSearch((ArrayList) data.getSample(), c, 0, pheIdx);
            calculate(modelName, pheIdx);
            count++;
        }
        bestSavedModelsMap.put(new Integer(or), savedModels);
    }

    public void mergeSearch(ArrayList subjects, String combination, int idxMarker, int pheIdx) {
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
                mergeSearch((ArrayList) subsets.get(key), com, idxMarker + 1, pheIdx);
            }
        } // if we've processed all attributes
        else {
            run(combination, subjects, pheIdx);
        }
    }

    protected void run(String com, ArrayList subsample, int pheIdx) {
        Suite s = new Suite(subsample, numTraits);
        s.summarize(pheIdx, offset[pheIdx]);
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

    public void calculate(String modelName, int pheIdx) {
        Set cellKeys = model.keySet();
        HeteroCombination hc = new HeteroCombination(numTraits, subdivision.getInterval());
            double thresholdRatio = 0;
            try {
                thresholdRatio = data.getTraitRatio(pheIdx, SNPIndex);
            } catch (DataFileException E) {
                E.printStackTrace(System.err);
            }
            for (int j = 0; j < subdivision.getInterval(); j++) {
                OneCVSet cvSet = new OneCVSet(j, modelName);
                Combination testingModels = (Combination) cvTestingSet.get(j);
                int status;
                Cell trCell;
                Cell tCell;
                for (Iterator e = cellKeys.iterator(); e.hasNext();) {
                    String cellKey = (String) e.next();
                    Suite fullSuite = (Suite) model.get(cellKey);
                    int fullposSubs = fullSuite.getPositiveSubjects(pheIdx);
                    int fullnegSubs = fullSuite.getNegativeSubjects(pheIdx);
                    double fullposScr = fullSuite.getPositiveScore(pheIdx);
                    double fullnegScr = fullSuite.getNegativeScore(pheIdx);
                    if (testingModels.containsKey(cellKey)) {
                        Suite testingSuite = (Suite) testingModels.get(cellKey);
                        int pos_Subs = testingSuite.getPositiveSubjects(pheIdx);
                        int neg_Subs = testingSuite.getNegativeSubjects(pheIdx);
                        double pos_Scr = testingSuite.getPositiveScore(pheIdx);
                        double neg_Scr = testingSuite.getNegativeScore(pheIdx);
                        status = ToolKit.Acertainment(fullposScr - pos_Scr, fullnegScr - neg_Scr, thresholdRatio);
                        trCell = new Cell(fullposSubs - pos_Subs, fullnegSubs - neg_Subs, fullposScr - pos_Scr, fullnegScr - neg_Scr, status);
                        tCell = new Cell(pos_Subs, neg_Subs, pos_Scr, neg_Scr, status);
                    } else {
                        status = ToolKit.Acertainment(fullposScr, fullnegScr, thresholdRatio);
                        trCell = new Cell(fullposSubs, fullnegSubs, fullposScr, fullnegScr, status);
                        tCell = new Cell(0, 0, 0, 0, -1);
                    }
                    cvSet.addTrainingModel(cellKey, trCell);
                    cvSet.addTestingModel(cellKey, tCell);
                }
                double trAccu = 0;
                try {
                    trAccu = ToolKit.Accuracy(cvSet.getTrainingSubdivision());
                } catch (ToolKitException E) {
                    E.printStackTrace(System.err);
                }
                cvSet.setStatistic(PublicData.TrainingAccuIdx, trAccu);
                double tAccu = 0;
                try {
                    tAccu = ToolKit.Accuracy(cvSet.getTestingSubdivision());
                } catch (ToolKitException E) {
                    E.printStackTrace(System.err);
                }
                cvSet.setStatistic(PublicData.TestingAccuIdx, tAccu);
                hc.SetStatistic(pheIdx, j, PublicData.TrainingAccuIdx, trAccu);
                hc.SetStatistic(pheIdx, j, PublicData.TestingAccuIdx, tAccu);
                if (count == 0) {
                    currBestStatistic[pheIdx][j] = trAccu;
                    bestKFoldResult.set(pheIdx, j, cvSet);
                    savedModels.save(modelName, model);
                } else if (currBestStatistic[pheIdx][j] < trAccu) {
                    currBestStatistic[pheIdx][j] = trAccu;
                    String keyAtIJ = bestKFoldResult.getKeyAt(pheIdx, j);
                    bestKFoldResult.set(pheIdx, j, cvSet);
                    savedModels.save(modelName, model);
                    savedModels.downsize(keyAtIJ);
                }
            }
        ModelHubs.put(modelName, hc);
    }

    public void summarise() {
        HashMap bestModelKeyMap = new HashMap();
        System.out.println("--------Order of interaction is " + bestKFoldResult.getOrder() + " --------");
        for (int i = 0; i < numTraits; i++) {
            System.out.println("========Statistics on trait " + (i + 1) + ".");
            for (int j = 0; j < subdivision.getInterval(); j++) {
                OneCVSet cvSet = bestKFoldResult.get(i, j);
                double testingAccu = 0;
                try {
                    testingAccu = ToolKit.Accuracy(cvSet.getTestingSubdivision());
                } catch (ToolKitException E) {
                    E.printStackTrace(System.err);
                }
                cvSet.setStatistic(PublicData.TestingAccuIdx, testingAccu);
            }
            double mean_TA = 0;
            double mean_TrA = 0;
            for (int j = 0; j < subdivision.getInterval(); j++) {
                OneCVSet cvSet = bestKFoldResult.get(i, j);
                double ta = cvSet.getStatistic(PublicData.TestingAccuIdx);
                double tra = cvSet.getStatistic(PublicData.TrainingAccuIdx);
                System.out.println("the chosen model at cv " + j + " is (" + cvSet.getCombination() + ")" + ", Testing accuracy is " + ta + ", Training accuracy is " + tra);
                mean_TA += ta;
                mean_TrA += tra;
            }
            mean_TA /= subdivision.getInterval();
            mean_TrA /= subdivision.getInterval();
            System.out.println("mean of Testing Accuracy is " + mean_TA);
            System.out.println("mean of Training Accuracy is " + mean_TrA);
            OneTraitResult i_thResult = (OneTraitResult) bestKFoldResult.get(i);
            i_thResult.summarise();
            System.out.println(i_thResult);
            bestModelKeyMap.put(new Integer(i), i_thResult.getBestModelKey());
        }
        SavedModels sModels = (SavedModels) bestSavedModelsMap.get(new Integer(order));
        sModels.saveModel(bestModelKeyMap);
    }

    public void print() {
        Set keys = (Set) ModelHubs.keySet();
        for (Iterator e = keys.iterator(); e.hasNext();) {
            String mn = (String) e.next();
            System.out.println(mn);
            HeteroCombination hc = (HeteroCombination) ModelHubs.get(mn);
            for (int i = 0; i < numTraits; i++) {
                for (int j = 0; j < subdivision.getInterval(); j++) {
                    System.out.println(hc.get(i, j, PublicData.TrainingAccuIdx) + "\t" + hc.get(i, j, PublicData.TestingAccuIdx));
                }
            }
        }
    }

    public void print(int or) {
        List com = comGenerator.get(or);
        for (Iterator e = com.iterator(); e.hasNext();) {
            String mn = (String) e.next();
            System.out.println(mn);
            HeteroCombination hc = (HeteroCombination) ModelHubs.get(mn);
            for (int i = 0; i < numTraits; i++) {
                for (int j = 0; j < subdivision.getInterval(); j++) {
                    System.out.println(hc.get(i, j, PublicData.TrainingAccuIdx) + "\t" + hc.get(i, j, PublicData.TestingAccuIdx));
                }
            }
        }
    }
}
