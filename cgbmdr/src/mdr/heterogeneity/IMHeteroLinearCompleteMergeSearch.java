package mdr.heterogeneity;

import java.io.*;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.HashMap;
import java.util.List;
import java.util.Set;

import algorithm.CombinationGenerator;
import algorithm.Subdivision;

import im.population.IMPopulation;
import im.GenomeScan;
import im.IMToolKit;
import im.IntervalPriorProbability;

import mdr.Cell;
import mdr.OneCVSet;

import publicAccess.PublicData;
import publicAccess.ToolKit;
import publicAccess.ToolKitException;

/**
 *
 * @author Guo-Bo Chen
 */
public class IMHeteroLinearCompleteMergeSearch extends AbstractHeteroMergeSearch {

    IMPopulation imp;
    GenomeScan gs;
    int[][] ChrInt;
    int[] MarkerIIPRowIdx;
    int currID;
    int phenoID;
    double peak;
    boolean hasEnvironment;
    ArrayList IDs;
    ArrayList permutatedIDs;
    ArrayList<HashMap<String, Double>> wholedataAccuracy;
    ArrayList<String> wholedataAccuracyKey;
    HashMap<String, IMSuite> trModel;
    HashMap<String, IMSuite> tModel;
    boolean searchSameChrInteraction;
    // indicating whether there's a chromosome containing more than one interval
    private boolean hasMoreThanOneInterval;

    public IMHeteroLinearCompleteMergeSearch(GenomeScan g, IMPopulation im, Subdivision sd, CombinationGenerator cg, int[] traits, double[] os, boolean ismooremdr) {
        super(sd, cg, traits.length, os, ismooremdr);
        gs = g;
        imp = im;
        IDs = im.getIDs();
        wholedataAccuracy = new ArrayList();

        for (int i = 0; i < numTraits; i++) {
            HashMap<String, Double> t = new HashMap();
            wholedataAccuracy.add(t);
        }
        peak = 0;
    }

    public void setSearch(boolean searchSameChromosome, boolean hasEnv) {
        searchSameChrInteraction = searchSameChromosome;
        hasEnvironment = hasEnv;
    }

    public void search(int or, int pheIdx, ArrayList sigComs) {
        if (testedInteraction.contains(new Integer(or))) {
            return;
        }
        wholedataAccuracyKey = new ArrayList();
        ChrInt = new int[or][2];
        MarkerIIPRowIdx = new int[or];
        order = or;

        for (Iterator e = sigComs.iterator(); e.hasNext();) {
            String modelName = (String) e.next();
            SNPIndex = ToolKit.StringToIntArray(modelName);
            ChrInt();
            if (!searchSameChrInteraction && goNext()) {
                continue;
            }
            int[] walks = new int[ChrInt.length];
            int[] currWalks = new int[ChrInt.length];
            for (int i = 0; i < ChrInt.length; i++) {
                IntervalPriorProbability iip = gs.getIPPTable(ChrInt[i][0], ChrInt[i][1]);
                walks[i] = iip.getWalks();
            }
            int[][] XCom = IMToolKit.XCombination(walks);
            int snpIdx = 0;
            for (int j = 0; j < XCom.length; j++) {
                StringBuffer tail = new StringBuffer(modelName);
                for (int k = 0; k < XCom[j].length; k++) {
                    tail.append(IMToolKit.separator + XCom[j][k]);
                }
                trModel = new HashMap();
                tModel = new HashMap();
                for (Iterator e1 = IDs.iterator(); e1.hasNext();) {
                    Integer id = (Integer) e1.next();
                    currID = id.intValue();
                    phenoID = (permutatedIDs == null) ? currID : ((Integer) (permutatedIDs.get(currID))).intValue();
                    for (int i = 0; i < SNPIndex.length; i++) {
                        MarkerIIPRowIdx[i] = gs.getIPPRowIndexForIndividual(currID, ChrInt[i][0], ChrInt[i][1]);
                    }
                    String combination = null;
                    if (hasEnvironment) {
                        combination = new String(imp.NestedEnvironment(currID));
                    } else {
                        combination = new String();
                    }
                    double p = 1;
                    mergeSearch(combination, 0, p, XCom[j], pheIdx);
                }
                wholedataAccuracyKey.add(tail.toString());
                calculateWholeDataAccuracy(tail.toString(), pheIdx);
            }
        }
    }

    public void search(int or, int pheIdx) {
        if (testedInteraction.contains(new Integer(or))) {
            return;
        }
        wholedataAccuracyKey = new ArrayList();
        ChrInt = new int[or][2];
        MarkerIIPRowIdx = new int[or];
        order = or;
        List com = (List) comGenerator.get(order);
        for (Iterator e = com.iterator(); e.hasNext();) {
            String modelName = (String) e.next();
            SNPIndex = ToolKit.StringToIntArray(modelName);
            ChrInt();
            if (!searchSameChrInteraction && goNext()) {
                continue;
            }
            int[] walks = new int[ChrInt.length];
            int[] currWalks = new int[ChrInt.length];
            for (int i = 0; i < ChrInt.length; i++) {
                IntervalPriorProbability iip = gs.getIPPTable(ChrInt[i][0], ChrInt[i][1]);
                walks[i] = iip.getWalks();
            }
            int[][] XCom = IMToolKit.XCombination(walks);
            int snpIdx = 0;
            for (int j = 0; j < XCom.length; j++) {
                StringBuffer tail = new StringBuffer(modelName);
                for (int k = 0; k < XCom[j].length; k++) {
                    tail.append(IMToolKit.separator + XCom[j][k]);
                }
                trModel = new HashMap();
                tModel = new HashMap();
                for (Iterator e1 = IDs.iterator(); e1.hasNext();) {
                    Integer id = (Integer) e1.next();
                    currID = id.intValue();
                    phenoID = (permutatedIDs == null) ? currID : ((Integer) (permutatedIDs.get(currID))).intValue();
                    for (int i = 0; i < SNPIndex.length; i++) {
                        MarkerIIPRowIdx[i] = gs.getIPPRowIndexForIndividual(currID, ChrInt[i][0], ChrInt[i][1]);
                    }
                    String combination = null;
                    if (hasEnvironment) {
                        combination = new String(imp.NestedEnvironment(currID));
                    } else {
                        combination = new String();
                    }
                    double p = 1;
                    mergeSearch(combination, 0, p, XCom[j], pheIdx);
                }
                wholedataAccuracyKey.add(tail.toString());
                calculateWholeDataAccuracy(tail.toString(), pheIdx);
            }
        }
    }

    private void mergeSearch(String combination, int idxMarker, double probability, int[] currWalks, int pheIdx) {
        if (idxMarker < SNPIndex.length) {
            IntervalPriorProbability iip = gs.getIPPTable(ChrInt[idxMarker][0], ChrInt[idxMarker][1]);
            ArrayList QTLtype = iip.QTLGenoType();
            for (int i = 0; i < QTLtype.size(); i++) {
                String qtl = (String) QTLtype.get(i);
                String com = new String();
                if (idxMarker == 0 && !hasEnvironment) {
                    com = qtl;
                } else {
                    com = combination + IMToolKit.separator + qtl;
                }
                mergeSearch(com, idxMarker + 1, probability * iip.PriorProbabilityAt(MarkerIIPRowIdx[idxMarker], currWalks[idxMarker], i), currWalks, pheIdx);
            }
        } else {
            run(combination, probability, pheIdx);
        }
    }

    private void run(String combination, double probability, int pheIdx) {
        IMSuite imsTr = trModel.get(combination);
        if (imsTr == null) {
            imsTr = new IMSuite(1);
            trModel.put(combination, imsTr);
        }
        imsTr.addScore(pheIdx, imp.ScoreAt(currID, pheIdx) * probability, probability);
        int intvl = ((Integer) dataPartitionMap.get(new Integer(currID))).intValue();
        IMSuite imsT = tModel.get(combination);
        if (imsT == null) {
            imsT = new IMSuite(1, subdivision.getInterval());
            tModel.put(combination, imsT);
        }
        imsT.addScore(pheIdx, intvl, imp.ScoreAt(currID, pheIdx) * probability, probability);
    }

    public void calculateWholeDataAccuracy(String modelName, int pheIdx) {
        Set cellKeys = trModel.keySet();
        double thresholdRatio = 0;
        thresholdRatio = imp.getTraitRatio(pheIdx);
        OneCVSet cvSet = new OneCVSet(0, modelName);
        for (Iterator e = cellKeys.iterator(); e.hasNext();) {
            String cellKey = (String) e.next();
            IMSuite fullSuite = (IMSuite) trModel.get(cellKey);
            int fullposSub = fullSuite.getCompletePositiveSubjects(pheIdx);
            int fullnegSub = fullSuite.getCompleteNegativeSubjects(pheIdx);
            double fullposScr = fullSuite.getCompletePositiveScore(pheIdx);
            double fullnegScr = fullSuite.getCompleteNegativeScore(pheIdx);
            int status = ToolKit.Acertainment(fullposScr, fullnegScr, thresholdRatio);
            Cell wholedataCell = new Cell(fullposSub, fullnegSub, fullposScr, fullnegScr, status);
            cvSet.addTrainingModel(cellKey, wholedataCell);
        }
        double wholedataAccu = 0;
        try {
            wholedataAccu = ToolKit.Accuracy(cvSet.getTrainingSubdivision());
        } catch (ToolKitException E) {
            E.printStackTrace(System.err);
        }
        HashMap<String, Double> t = (HashMap<String, Double>) wholedataAccuracy.get(pheIdx);
        t.put(modelName, new Double(wholedataAccu));
        if (wholedataAccu > peak) {
            peak = wholedataAccu;
        }
    }

    public void calculate(String modelName, int pheIdx) {
        Set cellKeys = trModel.keySet();
        HeteroCombination hc = new HeteroCombination(numTraits, subdivision.getInterval());
        hc.SetModelName(modelName);

        double thresholdRatio = 0;
        thresholdRatio = imp.getTraitRatio(pheIdx);
        for (int j = 0; j < subdivision.getInterval(); j++) {
            OneCVSet cvSet = new OneCVSet(j, modelName);
            int status;
            Cell trCell;
            Cell tCell;
            for (Iterator e = cellKeys.iterator(); e.hasNext();) {
                String cellKey = (String) e.next();
                IMSuite fullSuite = (IMSuite) trModel.get(cellKey);
                int fullposSubs = fullSuite.getCompletePositiveSubjects(pheIdx);
                int fullnegSubs = fullSuite.getCompleteNegativeSubjects(pheIdx);
                double fullposScr = fullSuite.getCompletePositiveScore(pheIdx);
                double fullnegScr = fullSuite.getCompleteNegativeScore(pheIdx);
                if (tModel.containsKey(cellKey)) {
                    IMSuite testingSuite = (IMSuite) tModel.get(cellKey);
                    int pos_Subs;
                    int neg_Subs;
                    double pos_Scr;
                    double neg_Scr;
                    if (testingSuite.ExistIt(pheIdx, j)) {
                        pos_Subs = testingSuite.getPositiveSubjects(pheIdx, j);
                        neg_Subs = testingSuite.getNegativeSubjects(pheIdx, j);
                        pos_Scr = testingSuite.getPositiveScore(pheIdx, j);
                        neg_Scr = testingSuite.getNegativeScore(pheIdx, j);
                        status = ToolKit.Acertainment(fullposScr - pos_Scr, fullnegScr - neg_Scr, thresholdRatio);
                        trCell = new Cell(fullposSubs - pos_Subs, fullnegSubs - neg_Subs, fullposScr - pos_Scr, fullnegScr - neg_Scr, status);
                        tCell = new Cell(pos_Subs, neg_Subs, pos_Scr, neg_Scr, status);
                    } else {
                        status = ToolKit.Acertainment(fullposScr, fullnegScr, thresholdRatio);
                        trCell = new Cell(fullposSubs, fullnegSubs, fullposScr, fullnegScr, status);
                        tCell = new Cell(0, 0, 0, 0, -1);
                    }
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
        }
        hc.summarise();
        ModelHubs.put(modelName, hc);
    }

    private void ChrInt() {
        int c = 0;
        int idx = 0;
        for (int i = 0; i < imp.ChromosomeNumber(); i++) {
            for (int j = 0; j < imp.IntervalNumberAtChromosome(i); j++) {
                if (c == SNPIndex[idx]) {
                    ChrInt[idx][0] = i;
                    ChrInt[idx][1] = j;
                    idx++;
                    if (idx == SNPIndex.length) {
                        break;
                    }
                }
                c++;
            }
            if (idx == SNPIndex.length) {
                break;
            }
        }
    }

    public void print() {
        System.out.println(wholedataAccuracyKey.size());
        HashMap<String, Double> t = (HashMap<String, Double>) wholedataAccuracy.get(0);

        for (int i = 0; i < wholedataAccuracyKey.size(); i++) {
            String key = wholedataAccuracyKey.get(i);
            if (!t.containsKey(key)) {
                System.out.println();
                continue;
            }
            String[] delim = key.split(IMToolKit.separator);
            StringBuffer m = new StringBuffer();
            for (int j = 0; j < order; j++) {
                m.append(delim[j]);
                if (j < (order - 1)) {
                    m.append(IMToolKit.separator);
                }
            }
            SNPIndex = ToolKit.StringToIntArray(m.toString());
            ChrInt();

            for (int j = 0; j < SNPIndex.length; j++) {
                double dis = imp.DistanceAt(ChrInt[j][0], ChrInt[j][1]);
                double dis1 = dis + Double.parseDouble(delim[j + order]) * gs.getStep();
                System.out.print(ChrInt[j][0] + "," + ChrInt[j][1] + "," + String.valueOf(dis1) + ",");
            }
            System.out.println(t.get(key));
        }
    }

    public void printBestEstimate(PrintStream out) {
        HashMap<String, Double> t = (HashMap<String, Double>) wholedataAccuracy.get(0);
        System.setOut(out);
        double bt = 0.0;
        String bKey = null;
        for (int i = 0; i < wholedataAccuracyKey.size(); i++) {
            String key = wholedataAccuracyKey.get(i);
            if (!t.containsKey(key)) {
                System.out.println();
                continue;
            }
            double b = ((Double) t.get(key)).doubleValue();

            String[] delim = key.split(IMToolKit.separator);
            StringBuffer m = new StringBuffer();
            for (int j = 0; j < order; j++) {
                m.append(delim[j]);
                if (j < (order - 1)) {
                    m.append(IMToolKit.separator);
                }
            }
            SNPIndex = ToolKit.StringToIntArray(m.toString());
            ChrInt();

            for (int j = 0; j < SNPIndex.length; j++) {
                double dis = imp.DistanceAt(ChrInt[j][0], ChrInt[j][1]);
                double dis1 = dis + Double.parseDouble(delim[j + order]) * gs.getStep();
                System.out.print(ChrInt[j][0] + " " + ChrInt[j][1] + " " + String.valueOf(dis1) + "\t");
            }
            System.out.println(t.get(key));

            if (b <= bt) {
                continue;
            } else {
                bt = b;
                bKey = key;
            }
        }

        String[] delim = bKey.split(IMToolKit.separator);
        StringBuffer m = new StringBuffer();
        for (int j = 0; j < order; j++) {
            m.append(delim[j]);
            if (j < (order - 1)) {
                m.append(IMToolKit.separator);
            }
        }
        SNPIndex = ToolKit.StringToIntArray(m.toString());
        ChrInt();
        System.out.print("Best estimation\t");
        for (int j = 0; j < SNPIndex.length; j++) {
            double dis = imp.DistanceAt(ChrInt[j][0], ChrInt[j][1]);
            double dis1 = dis + Double.parseDouble(delim[j + order]) * gs.getStep();
            System.out.print(ChrInt[j][0] + " " + ChrInt[j][1] + " " + String.valueOf(dis1) + "\t");
        }
        System.out.println(t.get(bKey));
        System.setOut(System.out);
    }

    private boolean goNext() {
        return hasMoreThanOneInterval;
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
