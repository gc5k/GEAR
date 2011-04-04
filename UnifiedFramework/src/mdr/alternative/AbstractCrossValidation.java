package mdr.alternative;

import java.util.AbstractList;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.TreeSet;

import mdr.Cell;
import mdr.alternative.DataSet;
import mdr.alternative.DataSetException;
import mdr.data.DataFile;
import algorithm.Partition;
import algorithm.Subdivision;
import mdr.Suite;

import mdr.Combination;

import publicAccess.ToolKit;
import publicAccess.ToolKitException;

/**
 *
 * @author Guo-Bo Chen
 */
public class AbstractCrossValidation extends AbstractList {

    protected ArrayList kFolderList = new ArrayList();

    protected boolean isMooreMDR;
    protected String cvKey;
    protected int[] scrIdx;
    protected int[] SNPIndex;
    protected double[] offset;

    protected ArrayList testingList;
    protected ArrayList testList;
    protected Combination model;
    protected HashMap fullModelResult;
    protected Set selectedPartition;

    protected DataSet data;
    protected Partition partitionInformation;
///    protected Subdivision subdivision;

    public AbstractCrossValidation(DataSet dr, Partition p, int[] idx, int[] sI, double[] os, Combination m, boolean ismoore) {
        data = dr;
///        subdivision = sub;
        partitionInformation = p;
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
        for (int i =0; i<partitionInformation.size(); i++) {
            Combination testingMap = new Combination();
            testingList.add(testingMap);
        }
/*
        testList = new ArrayList();
        for (int i=0; i<subdivision.getInterval(); i++) {
            Combination testMap = new Combination();
            testList.add(testMap);           
        }
*/
        for( int i=0; i<scrIdx.length; i++) {
            CVOneTraitResult k_thResult = new CVOneTraitResult();
            for( int j=0; j<testingList.size(); j++) {
///            for( int j=0; j<testList.size(); j++) {                
                CVPair cvPair = new CVPair();
                k_thResult.add(cvPair);
            }
            add(k_thResult);
        }      
    }

    public void testingCalculate(int sI, int division) {
        addSelectedPartition(division);
        CVOneTraitResult k_thTraitResult = (CVOneTraitResult) get(sI);
        CVPair cvPair = (CVPair) k_thTraitResult.get(division);
        double tAccu = 0;
        try {
            tAccu = ToolKit.Accuracy(cvPair.getTesting());
        } catch (ToolKitException E) {
            E.printStackTrace(System.err);
        }
        cvPair.setTestingAccuracy(tAccu);
    }

    public void calculate() {
        Set keys = model.keySet();
        for (Iterator e = keys.iterator(); e.hasNext();) {
            String key = (String) e.next();
            Suite fullSuite = (Suite) model.get(key);
            for (int i = 0; i < scrIdx.length; i++) {
                CVOneTraitResult k_thResult = (CVOneTraitResult) get(i);
                double thresholdRatio = 0;
                try {
                    thresholdRatio = data.getTraitRatio(SNPIndex, scrIdx[i]);
                } catch (DataSetException E) {
                    E.printStackTrace(System.err);
                }
                fullSuite.summarize(i, offset[i]);
                int fullModelStatus = ToolKit.Acertainment(fullSuite.getPositiveScore(i), fullSuite.getNegativeScore(i), thresholdRatio);
                Cell fullCell = new Cell(fullSuite.getPositiveSubjects(i), fullSuite.getNegativeSubjects(i),
                        fullSuite.getPositiveScore(i), fullSuite.getNegativeScore(i), fullModelStatus);
                fullModelResult.put(key, fullCell);
                int fullposSubs = fullSuite.getPositiveSubjects(i);
                int fullnegSubs = fullSuite.getNegativeSubjects(i);
                double fullposScr = fullSuite.getPositiveScore(i);
                double fullnegScr = fullSuite.getNegativeScore(i);
                for (int j = 0; j<testingList.size(); j++) {
                    CVPair cvPair = (CVPair) k_thResult.get(j);
                    Combination testingModels = (Combination) testingList.get(j);
                    int status;
                    Cell trCell;
                    Cell tCell;
                    if (testingModels.containsKey(key)) {
                        Suite testingSuite = (Suite) testingModels.get(key);
                        testingSuite.summarize(i, offset[i]);
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
/*
    public void kFolderSubdivision() {
        Set keys = model.keySet();
        HashMap division = subdivision.getDivision();
        for (Iterator e = keys.iterator(); e.hasNext();) {
            String key = (String) e.next();
            Suite fullSuite = (Suite) model.get(key);
            for (Iterator e1 = fullSuite.iterator(); e1.hasNext();) {
                DataSet.Subject sub = (DataSet.Subject) e1.next();
                Integer ID = sub.getIntegerID();
                int d = ((Integer) division.get(ID)).intValue();
                Combination suiteMap = (Combination) testList.get(d);
                Suite S = (Suite) suiteMap.get(key);
                if (S == null) {
                    S = new Suite();
                    suiteMap.put(key, S);
                }
                S.add(sub);
            }
            
            for (int i = 0; i < scrIdx.length; i++) {
                CVOneTraitResult k_thTraitResult = (CVOneTraitResult) get(i);
                double thresholdRatio = 0;
                try {
                    thresholdRatio = data.getTraitRatio(SNPIndex, scrIdx[i]);
                } catch (DataSetException E) {
                    E.printStackTrace(System.err);
                }
                fullSuite.summarize(i, offset[i]);

                int fullModelStatus = ToolKit.acertainment(fullSuite.getPositiveScore(), fullSuite.getNegativeScore(), thresholdRatio);
                Cell fullCell = new Cell(fullSuite.getPositiveSubjects(), fullSuite.getNegativeSubjects(),
                        fullSuite.getPositiveScore(), fullSuite.getNegativeScore(), fullModelStatus);
                fullModelResult.put(key, fullCell);
                int fullposSubs = fullSuite.getPositiveSubjects();
                int fullnegSubs = fullSuite.getNegativeSubjects();
                double fullposScr = fullSuite.getPositiveScore();
                double fullnegScr = fullSuite.getNegativeScore();
                for (int j = 0; j<testList.size(); j++) {
                    CVPair cvPair = (CVPair) k_thTraitResult.get(j);
                    Combination testingModel = (Combination) testList.get(j);
                    int status;
                    Cell trCell;
                    Cell tCell;
                    if (testingModel.containsKey(key)) {
                        Suite testingSuite = (Suite) testingModel.get(key);
                        testingSuite.summarize(scrIdx[i], offset[i]);
                        int pos_Subs = testingSuite.getPositiveSubjects();
                        int neg_Subs = testingSuite.getNegativeSubjects();
                        double pos_Scr = testingSuite.getPositiveScore();
                        double neg_Scr = testingSuite.getNegativeScore();
                        status = ToolKit.acertainment(fullposScr - pos_Scr, fullnegScr - neg_Scr, thresholdRatio);
                        trCell = new Cell(fullposSubs - pos_Subs, fullnegSubs - neg_Subs, fullposScr - pos_Scr, fullnegScr - neg_Scr, status);
                        tCell = new Cell(pos_Subs, neg_Subs, pos_Scr, neg_Scr, status);
                    } else {
                        status = ToolKit.acertainment(fullposScr, fullnegScr, thresholdRatio);
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
*/
    public void kFolder() {
        Set keys = model.keySet();
        for (Iterator e = keys.iterator(); e.hasNext();) {
            String key = (String) e.next();
            Suite fullSuite = (Suite) model.get(key);
            ArrayList temp_fullSuite = fullSuite.getSubjects();
            for (Iterator e1 = temp_fullSuite.iterator(); e1.hasNext();) {
                DataFile.Subject sub = (DataFile.Subject) e1.next();
                Integer ID = sub.getIntegerID();
                int count = 0;
                for (Iterator e2 = partitionInformation.iterator(); e2.hasNext();) {
                    HashSet IDs = (HashSet) e2.next();
                    if (IDs.contains(ID)) {
                        Combination suiteMap = (Combination) testingList.get(count);
                        Suite S = (Suite) suiteMap.get(key);
                        if (S == null) {
                            S = new Suite(scrIdx.length);
                            suiteMap.put(key, S);
                        }
                        S.add(sub);
                    }
                    count++;
                }
            }
        }
    }

    public void printKFolderData() {
        int c = 0;
        for (Iterator e = testingList.iterator(); e.hasNext();) {
//        for (Iterator e = testList.iterator(); e.hasNext();) {        
            System.out.println("C=" + c++);
            HashMap tstMap = (HashMap) e.next();
            TreeSet tstkeyset = new TreeSet(tstMap.keySet());
            for (Iterator se = tstkeyset.iterator(); se.hasNext();) {
                String key = (String) se.next();
                System.out.println("testing " + key);
                Suite tstSuite = (Suite) tstMap.get(key);
                ArrayList temp_tstSuite = tstSuite.getSubjects();
                for (Iterator le = temp_tstSuite.iterator(); le.hasNext();) {
                    DataSet.Subject sub = (DataSet.Subject) le.next();
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

    public class CVOneTraitResult extends AbstractList {

        private ArrayList ResultList = new ArrayList();

        public CVOneTraitResult() {
        }

        public void add(int idx, Object o) {
            modCount++;
            ResultList.add(idx, o);
        }

        public Object get(int idx) {
            return ResultList.get(idx);
        }

        public void print() {
            int cv = 0;
            for (Iterator e = ResultList.iterator(); e.hasNext();) {
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

        public Object remove(int idx) {
            modCount++;
            return ResultList.remove(idx);
        }

        public int size() {
            return ResultList.size();
        }
    }    
}