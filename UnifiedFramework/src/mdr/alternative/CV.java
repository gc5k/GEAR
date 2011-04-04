package mdr.alternative;

import mdr.Combination;
import mdr.*;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.TreeSet;

import mdr.Cell;
import mdr.data.DataFile;
import mdr.DataFileException;
import algorithm.Subdivision;

import mdr.Suite;
import publicAccess.ToolKit;
import publicAccess.ToolKitException;

/**
 *
 * @author Guo-Bo Chen
 */
public class CV extends AbstractList {

    protected ArrayList kFolderList = new ArrayList();

    protected boolean isMooreMDR;
    protected String cvKey;
    protected int[] scrIdx;
    protected int[] SNPIndex;
    protected double[] offset;

    protected ArrayList testingList;
    protected Combination model;
    protected HashMap fullModelResult;
    protected Set selectedPartition;

    protected DataFile data;
    protected Subdivision subdivision;

    public void calculate() {
        Set keys = model.keySet();
        for (Iterator e = keys.iterator(); e.hasNext();) {
            String key = (String) e.next();
            Suite fullSuite = (Suite) model.get(key);
            for (int i = 0; i < scrIdx.length; i++) {
                CVOneTraitResult k_thResult = (CVOneTraitResult) get(i);
                double thresholdRatio = 0;
                try {
                    thresholdRatio = data.getTraitRatio(i, SNPIndex);
                } catch (DataFileException E) {
                    E.printStackTrace(System.err);
                }
                fullSuite.summarize(scrIdx[i], offset[i]);
                int fullposSubs = fullSuite.getPositiveSubjects(i);
                int fullnegSubs = fullSuite.getNegativeSubjects(i);
                double fullposScr = fullSuite.getPositiveScore(i);
                double fullnegScr = fullSuite.getNegativeScore(i);
                int fullModelStatus = ToolKit.Acertainment(fullposScr, fullnegScr, thresholdRatio);
                Cell fullCell = new Cell(fullposSubs, fullnegSubs,
                        fullposScr, fullnegScr, fullModelStatus);
                fullModelResult.put(key, fullCell);
                for (int j = 0; j<testingList.size(); j++) {
                    CVPair cvPair = (CVPair) k_thResult.get(j);
                    Combination testingModels = (Combination) testingList.get(j);
                    int status;
                    Cell trCell;
                    Cell tCell;
                    if (testingModels.containsKey(key)) {
                        Suite testingSuite = (Suite) testingModels.get(key);
                        testingSuite.summarize(scrIdx[i], offset[i]);
                        int pos_Subs = testingSuite.getPositiveSubjects(i);
                        int neg_Subs = testingSuite.getNegativeSubjects(i);
                        double pos_Scr = testingSuite.getPositiveScore(i);
                        double neg_Scr = testingSuite.getNegativeScore(i);
                        status = ToolKit.Acertainment(fullposScr - pos_Scr, fullnegScr - neg_Scr, thresholdRatio);
                        trCell = new Cell(fullposSubs - pos_Subs, fullnegSubs - neg_Subs, fullposScr - pos_Scr, fullnegScr - neg_Scr, status);
                        tCell = new Cell(pos_Subs, neg_Subs, pos_Scr, neg_Scr, status);
                    } else {
                        status = ToolKit.Acertainment(fullposScr, fullnegScr, thresholdRatio);
                        trCell = new Cell(fullposSubs, fullnegSubs, fullposScr, fullnegScr, status);
                        tCell = new Cell(0, 0, 0, 0, -1);
                    }
                    cvPair.addTrainingModel(key, trCell);
                    cvPair.addTestingModel(key, tCell);
                }
            }
        }

        for (int i = 0; i < scrIdx.length; i++) {
            CVOneTraitResult k_thResult = (CVOneTraitResult) get(i);
            for (int j = 0; j < testingList.size(); j++) {
                CVPair cvPair = (CVPair) k_thResult.get(j);
                double trAccu = 0;
                try {
                    trAccu = ToolKit.Accuracy(cvPair.getTraining());
                } catch (ToolKitException E) {
                    E.printStackTrace(System.err);
                }
                cvPair.setTrainingAccuracy(trAccu);
                if (!isMooreMDR) {
                    testingCalculate(i, j);
                }
            }
        }
    }

    public void testingCalculate(int sI, int subdivision) {
        addSelectedPartition(subdivision);
        CVOneTraitResult k_thResult = (CVOneTraitResult) get(sI);
        CVPair cvPair = (CVPair) k_thResult.get(subdivision);
        double tAccu = 0;
        try {
            tAccu = ToolKit.Accuracy(cvPair.getTesting());
        } catch (ToolKitException E) {
            E.printStackTrace(System.err);
        }
        cvPair.setTestingAccuracy(tAccu);
    }

    public CV(DataFile dr, Subdivision sd, int[] idx, int[] sI, double[] os, Combination m, boolean ismoore) {
        data = dr;
        subdivision = sd;
        SNPIndex = new int[idx.length];
        System.arraycopy(idx, 0, SNPIndex, 0, idx.length);
        scrIdx = new int[sI.length];
        System.arraycopy(sI, 0, scrIdx, 0, sI.length);
        cvKey = ToolKit.IntArrayToString(SNPIndex);
        offset = new double[os.length];
        System.arraycopy(os, 0, offset, 0, os.length);
        model = m;
        selectedPartition = new HashSet();
        fullModelResult = new HashMap();
        isMooreMDR = ismoore;
        testingList = new ArrayList();
        for (int i =0; i<subdivision.getInterval(); i++) {
            Combination testingMap = new Combination();
            testingList.add(testingMap);
        }

        for( int i=0; i<scrIdx.length; i++) {
            CVOneTraitResult k_thResult = new CVOneTraitResult();
            for( int j=0; j<testingList.size(); j++) {
                CVPair cvPair = new CVPair();
                k_thResult.add(cvPair);
            }
            add(k_thResult);
        }        
    }

    public void kFolder() {
        Set keys = model.keySet();
        HashMap divisionMap = subdivision.getDivision();
        for (Iterator e = keys.iterator(); e.hasNext();) {
            String key = (String) e.next();
            Suite fullSuite = (Suite) model.get(key);
            ArrayList temp_fullSuite = fullSuite.getSubjects();
            for (Iterator e1 = temp_fullSuite.iterator(); e1.hasNext();) {
                DataFile.Subject sub = (DataFile.Subject) e1.next();
                Integer ID = sub.getIntegerID();
                int d = ((Integer) divisionMap.get(ID)).intValue();
                Combination suiteMap = (Combination) testingList.get(d);
                Suite S = (Suite) suiteMap.get(key);
                if (S == null) {
                    S = new Suite(scrIdx.length);
                    suiteMap.put(key, S);
                }
                S.add(sub);
            }
        }
    }


    public void kFolderSubdivision() {
        Set keys = model.keySet();
        HashMap division = subdivision.getDivision();
        for (Iterator e = keys.iterator(); e.hasNext();) {
            String key = (String) e.next();
            Suite fullSuite = (Suite) model.get(key);
            ArrayList temp_fullSuite = fullSuite.getSubjects();
            for (Iterator e1 = temp_fullSuite.iterator(); e1.hasNext();) {
                DataFile.Subject sub = (DataFile.Subject) e1.next();
                Integer ID = sub.getIntegerID();
                int d = ((Integer) division.get(ID)).intValue();
                Combination suiteMap = (Combination) testingList.get(d);
                Suite S = (Suite) suiteMap.get(key);
                if (S == null) {
                    S = new Suite(scrIdx.length);
                    suiteMap.put(key, S);
                }
                S.add(sub);
            }
            
            for (int i = 0; i < scrIdx.length; i++) {
                CVOneTraitResult k_thTraitResult = (CVOneTraitResult) get(i);
                double thresholdRatio = 0;
                try {
                    thresholdRatio = data.getTraitRatio(i, SNPIndex);
                } catch (DataFileException E) {
                    E.printStackTrace(System.err);
                }
                fullSuite.summarize(scrIdx[i], offset[i]);
                int fullposSubs = fullSuite.getPositiveSubjects(i);
                int fullnegSubs = fullSuite.getNegativeSubjects(i);
                double fullposScr = fullSuite.getPositiveScore(i);
                double fullnegScr = fullSuite.getNegativeScore(i);
                int fullModelStatus = ToolKit.Acertainment(fullposScr, fullnegScr, thresholdRatio);
                Cell fullCell = new Cell(fullposSubs, fullnegSubs,
                        fullposScr, fullnegScr, fullModelStatus);
                fullModelResult.put(key, fullCell);

                for (int j = 0; j<testingList.size(); j++) {
                    CVPair cvPair = (CVPair) k_thTraitResult.get(j);
                    Combination testingModel = (Combination) testingList.get(j);
                    int status;
                    Cell trCell;
                    Cell tCell;
                    if (testingModel.containsKey(key)) {
                        Suite testingSuite = (Suite) testingModel.get(key);
                        testingSuite.summarize(scrIdx[i], offset[i]);
                        int pos_Subs = testingSuite.getPositiveSubjects(i);
                        int neg_Subs = testingSuite.getNegativeSubjects(i);
                        double pos_Scr = testingSuite.getPositiveScore(i);
                        double neg_Scr = testingSuite.getNegativeScore(i);
                        status = ToolKit.Acertainment(fullposScr - pos_Scr, fullnegScr - neg_Scr, thresholdRatio);
                        trCell = new Cell(fullposSubs - pos_Subs, fullnegSubs - neg_Subs, fullposScr - pos_Scr, fullnegScr - neg_Scr, status);
                        tCell = new Cell(pos_Subs, neg_Subs, pos_Scr, neg_Scr, status);
                    } else {
                        status = ToolKit.Acertainment(fullposScr, fullnegScr, thresholdRatio);
                        trCell = new Cell(fullposSubs, fullnegSubs, fullposScr, fullnegScr, status);
                        tCell = new Cell(0, 0, 0, 0, -1);
                    }
                    cvPair.addTrainingModel(key, trCell);
                    cvPair.addTestingModel(key, tCell);
                    double trAccu = 0;
                    try {
                        trAccu = ToolKit.Accuracy(cvPair.getTraining());
                    } catch (ToolKitException E) {
                        E.printStackTrace(System.err);
                    }
                    cvPair.setTrainingAccuracy(trAccu);
                    if (!isMooreMDR) {
                        testingCalculate(i, j);
                    }
                }
            }
        }
    }
    
    
    public void printKFolderData() {
        int c = 0;
        for (Iterator e = testingList.iterator(); e.hasNext();) {
            System.out.println("C=" + c++);
            HashMap tstMap = (HashMap) e.next();
            TreeSet tstkeyset = new TreeSet(tstMap.keySet());
            for (Iterator se = tstkeyset.iterator(); se.hasNext();) {
                String key = (String) se.next();
                System.out.println("testing " + key);
                Suite tstSuite = (Suite) tstMap.get(key);
                ArrayList temp_tstSuite = tstSuite.getSubjects();
                for (Iterator le = temp_tstSuite.iterator(); le.hasNext();) {
                    DataFile.Subject sub = (DataFile.Subject) le.next();
                    System.out.println(sub);
                }
            }
        }
    }

    public void add(int idx, Object o) {
        modCount++;
        kFolderList.add(idx, o);
    }

    public Object get(int idx) {
        return kFolderList.get(idx);
    }

    public String getCVKey() {
        return cvKey;
    }
    
    public void print() {
        printKFolderData();
    }

    public int size() {
        return kFolderList.size();
    }
    
    protected void addSelectedPartition(int i) {
        selectedPartition.add(new Integer(i));
    }
    
    public void testPrint() {
        for(Iterator e = kFolderList.iterator(); e.hasNext(); ) {
            CVOneTraitResult cvResult = (CVOneTraitResult) e.next();
            cvResult.print();
        }
    }

    public boolean isMooreMDR() {
        return isMooreMDR;
    }    

    public class CVPair {

        HashMap trainingMap;
        HashMap testingMap;
        double testingAccuracy;
        double trainingAccuracy;
        
        public CVPair(HashMap trn, HashMap tst) {
            trainingMap = trn;
            testingMap = tst;
        }

        public double getTestingAccuracy() {
            return testingAccuracy;
        }
        
        public void setTestingAccuracy(double tstAccu) {
            testingAccuracy = tstAccu;
        }
        
        public double getTrainingAccuracy() {
            return trainingAccuracy;
        }
        public void setTrainingAccuracy(double trAccu) {
            trainingAccuracy = trAccu;
        }     

        public CVPair() {
            trainingMap = new HashMap();
            testingMap = new HashMap();
        }
        
        public void addTrainingModel(String key, Object o) {
            trainingMap.put(key, o);
        }
        
        public void addTestingModel(String key, Object o) {
            testingMap.put(key, o);
        }
        
        public HashMap getTraining() {
            return trainingMap;
        }

        public HashMap getTesting() {
            return testingMap;
        }
    }
    
    public class CVOneTraitResult extends ArrayList {

        public void print() {
            int cv = 0;
            for (Iterator e = this.iterator(); e.hasNext();) {
                CVPair cvPair = (CVPair) e.next();
                System.out.println("CV result " + cv + ", Tesint Accuracy " + cvPair.getTestingAccuracy() + ", Training Accuracy " + cvPair.getTrainingAccuracy());
                HashMap trMap = cvPair.getTraining();
                TreeSet trkeyset = new TreeSet(trMap.keySet());
                for (Iterator se = trkeyset.iterator(); se.hasNext();) {
                    String key = (String) se.next();
//                    System.out.println(key);
                    Cell trCell = (Cell) trMap.get(key);
//                    System.out.println(trCell);
                    if (!isMooreMDR || selectedPartition.contains(new Integer(cv))) {
                        HashMap tMap = cvPair.getTesting();
                        Cell tCell = (Cell) tMap.get(key);
//                            System.out.println(tCell);
                    }
                }
                cv++;
            }
        }
    }    
}