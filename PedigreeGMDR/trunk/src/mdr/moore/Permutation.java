package mdr.moore;

import mdr.Combination;
import mdr.Cell;
import mdr.Suite;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.HashMap;
import java.util.Random;
import java.util.Set;

import mdr.data.DataFile;
import mdr.DataFileException;
import algorithm.Subdivision;

import mdr.SavedModels;
import mdr.OneCVSet;

import publicAccess.PublicData;
import publicAccess.ToolKit;
import publicAccess.ToolKitException;

/**
 *
 * @author Guo-Bo Chen
 */
public class Permutation {

    private int numTraits;
    private double[] offset;
    private int[] SNPIndex;
    private String[] modelKey;
    private int replication;

    private ArrayList Trait = new ArrayList();
    private ArrayList TraitIndex = new ArrayList();    
    private HashMap dataPartitionMap;
    private ArrayList<PermutationResult> permutationResult = new ArrayList();
    private ArrayList<PCombination> cvTestingSet = new ArrayList();//always keeps the K testing models of the current combination 
    private ArrayList<PCombination> pcomList = new ArrayList();
    
    private DataFile data;
    private SavedModels savedModels;
    private Subdivision subdivision;

    public Permutation(DataFile dr, SavedModels sm, double[] os, Subdivision sub, int rep) {

        data = dr;
        numTraits = sm.size();
        offset = new double[os.length];
        System.arraycopy(os, 0, offset, 0, os.length);

        modelKey = new String[numTraits];
        savedModels = sm;
        subdivision = sub;
        replication = rep + 1;
        for( int i = 0; i < numTraits; i++) {
            modelKey[i] = sm.getModelKey(new Integer(i));
            Trait.add(new ArrayList());
            TraitIndex.add(new HashMap());
            permutationResult.add(new PermutationResult(replication, subdivision.getInterval()));
            pcomList.add(new PCombination(numTraits));
        }

        Combination com = (Combination) savedModels.getBestModel(new Integer(0));
        dataPartitionMap = sub.getDivision();

        for (int i = 0; i < subdivision.getInterval(); i++) {
            PCombination testingMap = new PCombination(numTraits);
            cvTestingSet.add(testingMap);
        }

        initial(sm);
    }

    private void initial(SavedModels sm) {
        for( int i = 0; i < sm.size(); i++) {
            Combination com = sm.getBestModel(new Integer(i));
            Set keys = com.keySet();
            PCombination pcombination = (PCombination) pcomList.get(i);
            for (Iterator e = keys.iterator(); e.hasNext();) {
                String key = (String) e.next();
                Suite s = (Suite) com.get(key);
                ArrayList subjects = s.getSubjects();
                ArrayList trait = (ArrayList) Trait.get(i);
                HashMap traitIndex = (HashMap) TraitIndex.get(i);
                for (Iterator e1 = subjects.iterator(); e1.hasNext();) {
                    DataFile.Subject sub = (DataFile.Subject) e1.next();
                    Double subscore = sub.getDoubleScore(i);
                    if (subscore != null) {
                        trait.add(subscore);
                        traitIndex.put(sub.getIntegerID(), new Integer(trait.size() - 1));
                    }
                }
                PSuite ps = new PSuite(numTraits);
                for (Iterator e1 = subjects.iterator(); e1.hasNext();) {
                    DataFile.Subject sub = (DataFile.Subject) e1.next();
                    ps.addSubject(sub, TraitIndex);
                }
                pcombination.addPSuite(key, ps);
            }
        }
        calculateFullModel();        
    }

    private void calculateFullModel() {
        for (int i = 0; i < pcomList.size(); i++) {
            PCombination pcombination = pcomList.get(i);
            Set keys = pcombination.keySet();
            for (Iterator e1 = keys.iterator(); e1.hasNext();) {
                String key = (String) e1.next();
                PSuite ps = (PSuite) pcombination.get(key);
                ps.summarize(i, offset[i], Trait);
            }
        }
    }

    public void evaluate() {
        assignKFold();
        calculate();
    }

    public PermutationResult getResult(int idx) {
        return (PermutationResult) permutationResult.get(idx);
    }

    protected void assignKFold() {
        for (int i = 0; i < pcomList.size(); i++) {
            PCombination pcombination = pcomList.get(i);
            Set keys = pcombination.keySet();
            for (Iterator e = keys.iterator(); e.hasNext();) {
                String key = (String) e.next();
                ArrayList subsample = ((PSuite) pcombination.get(key)).getSubjects();
                for (Iterator e1 = subsample.iterator(); e1.hasNext();) {
                    DataFile.Subject sub = (DataFile.Subject) e1.next();
                    Integer ID = sub.getIntegerID();
                    int d = ((Integer) dataPartitionMap.get(ID)).intValue();
                    PCombination pcom = (PCombination) cvTestingSet.get(d);
                    PSuite PS = (PSuite) pcom.get(key);
                    if (PS == null) {
                        PS = new PSuite(numTraits);
                        pcom.put(key, PS);
                    }
                    PS.addSubject(sub, TraitIndex);
                }
            }
        }
    }

    public void calculate() {
        for (int k = 0; k < replication; k++) {
            for (int i = 0; i < numTraits; i++) {
                PCombination pcombination = pcomList.get(i);
                Set cellKeys = pcombination.keySet();
                if ( k > 0) {
                    shuffle(i, k);
                }
                calculateFullModel();
                double thresholdRatio = 0;
                SNPIndex = ToolKit.StringToIntArray(modelKey[i]);
                try {
                    thresholdRatio = data.getTraitRatio(i, SNPIndex);
                } catch (DataFileException E) {
                    E.printStackTrace(System.err);
                }
                for (int j = 0; j < subdivision.getInterval(); j++) {
                    OneCVSet cvSet = new OneCVSet(j, modelKey[i]);
                    PCombination testingModels = (PCombination) cvTestingSet.get(j);
                    int status;
                    Cell trCell;
                    Cell tCell;
                    for (Iterator e = cellKeys.iterator(); e.hasNext();) {
                        String cellKey = (String) e.next();
                        PSuite fullPSuite = (PSuite) pcombination.get(cellKey);
                        int fullposSubs = fullPSuite.getPositiveSubjects(i);
                        int fullnegSubs = fullPSuite.getNegativeSubjects(i);
                        double fullposScr = fullPSuite.getPositiveScore(i);
                        double fullnegScr = fullPSuite.getNegativeScore(i);
                        if (testingModels.containsKey(cellKey)) {
                            PSuite testingPSuite = (PSuite) testingModels.get(cellKey);
                            testingPSuite.summarize(i, offset[i], Trait);
                            int pos_Subs = testingPSuite.getPositiveSubjects(i);
                            int neg_Subs = testingPSuite.getNegativeSubjects(i);
                            double pos_Scr = testingPSuite.getPositiveScore(i);
                            double neg_Scr = testingPSuite.getNegativeScore(i);
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
                    ((PermutationResult) (permutationResult.get(i))).set(k, PublicData.TrainingAccuIdx, j, trAccu);
                    double tAccu = 0;
                    try {
                        tAccu = ToolKit.Accuracy(cvSet.getTestingSubdivision());
                    } catch (ToolKitException E) {
                        E.printStackTrace(System.err);
                    }
                    ((PermutationResult) (permutationResult.get(i))).set(k, PublicData.TestingAccuIdx, j, tAccu);
                }
            }
        }
    }

    public ArrayList getPCombination() {
        return pcomList;
    }

    protected void shuffle(int idx, long seed) {
        ArrayList trait = (ArrayList) Trait.get(idx);
        Random rnd = new Random(seed);
        Collections.shuffle(trait, rnd);
    }
}
