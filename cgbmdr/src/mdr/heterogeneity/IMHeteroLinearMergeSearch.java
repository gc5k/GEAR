package mdr.heterogeneity;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.HashMap;
import java.util.HashSet;
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
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class IMHeteroLinearMergeSearch extends AbstractHeteroMergeSearch {

    IMPopulation imp;
    GenomeScan gs;
    int[][] ChrInt;
    int[] MarkerIIPRowIdx;
    double[] bestTestingAccu;
    double peak;
    int currID;
    int phenoID;
    boolean hasEnvironment;
    ArrayList IDs;
    ArrayList permutatedIDs;
    ArrayList<HashMap<String, Double>> wholedataAccuracy;
    HashMap<String, IMSuite> trModel;
    HashMap<String, IMSuite> tModel;
    HashMap<Integer, String> bestModel;
    boolean searchSameChrInteraction;
    // indicating whether there's a chromosome containing more than one interval
    private boolean hasMoreThanOneInterval;

    public IMHeteroLinearMergeSearch(GenomeScan g, IMPopulation im, Subdivision sd, CombinationGenerator cg, int[] traits, double[] os, boolean ismooremdr) {
        super(sd, cg, traits.length, os, ismooremdr);
        gs = g;
        imp = im;
        IDs = im.getIDs();
        searchSameChrInteraction = true;
        wholedataAccuracy = new ArrayList();
        for (int i = 0; i < numTraits; i++) {
            HashMap<String, Double> t = new HashMap();
            wholedataAccuracy.add(t);
        }
        peak = 0;
    }

    public IMHeteroLinearMergeSearch(GenomeScan g, IMPopulation im, Subdivision sd, CombinationGenerator cg, int[] traits, double[] os, boolean ismooremdr, ArrayList ids) {
        super(sd, cg, traits.length, os, ismooremdr);
        gs = g;
        imp = im;
        IDs = im.getIDs();
        permutatedIDs = ids;
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

    /**
     *
     * @param or    interaction order
     */
    public void search(int or, int pheIdx) {
        if (testedInteraction.contains(new Integer(or))) {
            return;
        }
        ChrInt = new int[or][2];
        MarkerIIPRowIdx = new int[or];
        order = or;
        List com = (List) comGenerator.get(order);
        bestModel = new HashMap();
        bestTestingAccu = new double[subdivision.getInterval()];
        for (Iterator e = com.iterator(); e.hasNext();) {
            String modelName = (String) e.next();
            SNPIndex = ToolKit.StringToIntArray(modelName);
            ChrInt();
            if (!searchSameChrInteraction && goNext()) {
                continue;
            }
            trModel = new HashMap();
            tModel = new HashMap();
            HashMap realizedCom = new HashMap();
            for (Iterator e1 = IDs.iterator(); e1.hasNext();) {
                Integer id = (Integer) e1.next();
                currID = id.intValue();
                phenoID = (permutatedIDs == null) ? currID : ((Integer) (permutatedIDs.get(currID))).intValue();
                StringBuffer sb = new StringBuffer();
                String combination = null;
                if (hasEnvironment) {
                    combination = new String(imp.NestedEnvironment(currID));
                    sb.append(imp.NestedEnvironment(currID));
                } else {
                    combination = new String();
                }
                for (int i = 0; i < SNPIndex.length; i++) {
                    MarkerIIPRowIdx[i] = gs.getIPPRowIndexForIndividual(currID, ChrInt[i][0], ChrInt[i][1]);
                    sb.append(imp.MarkerAt(currID, ChrInt[i][0], ChrInt[i][1]) + IMToolKit.separator + imp.MarkerAt(currID, ChrInt[i][0], ChrInt[i][1]) + IMToolKit.separator);
                }
                double p = 1;
                mergeSearch(combination, 0, p, pheIdx);
            }
            calculateWholeDataAccuracy(modelName, pheIdx);
        }
    }

    public void search(int or, int pheIdx, int[][] ci) {
        if (testedInteraction.contains(new Integer(or))) {
            return;
        }
        ChrInt = new int[or][2];
        MarkerIIPRowIdx = new int[or];
        order = or;
        List com = (List) comGenerator.get(order);
        bestModel = new HashMap();
        bestTestingAccu = new double[subdivision.getInterval()];
        for (Iterator e = com.iterator(); e.hasNext();) {
            String modelName = (String) e.next();
            SNPIndex = ToolKit.StringToIntArray(modelName);
            ChrInt();
            boolean flag = true;
            for (int i = 0; i < ci.length; i++) {
                if (ci[i][0] != ChrInt[i][0] || ci[i][1] != ChrInt[i][1]) {
                    flag = false;
                    break;
                }
            }
            if (!flag) {
                continue;
            }
            if (!searchSameChrInteraction && goNext()) {
                continue;
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
                mergeSearch(combination, 0, p, pheIdx);
            }
            calculateWholeDataAccuracy(modelName, pheIdx);
//            calculate(modelName);
        }
    }

    private void mergeSearch(String combination, int idxMarker, double probability, int pheIdx) {
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
                mergeSearch(com, idxMarker + 1, probability * iip.MidPriorProbability(MarkerIIPRowIdx[idxMarker], i), pheIdx);
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
        double scr = imp.ScoreAt(phenoID, pheIdx);
        imsTr.addScore(pheIdx, scr * probability, probability);
        int intvl = ((Integer) dataPartitionMap.get(new Integer(currID))).intValue();
        IMSuite imsT = tModel.get(combination);
        if (imsT == null) {
            imsT = new IMSuite(1, subdivision.getInterval());
            tModel.put(combination, imsT);
        }
        imsT.addScore(pheIdx, intvl, scr * probability, probability);
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

    public void calculate(String modelName) {
        Set cellKeys = trModel.keySet();
        HeteroCombination hc = new HeteroCombination(numTraits, subdivision.getInterval());
        hc.SetModelName(modelName);
        for (int i = 0; i < numTraits; i++) {
            double thresholdRatio = 0;
            thresholdRatio = imp.getTraitRatio(i);
            for (int j = 0; j < subdivision.getInterval(); j++) {
                OneCVSet cvSet = new OneCVSet(j, modelName);
                int status;
                Cell trCell;
                Cell tCell;
                for (Iterator e = cellKeys.iterator(); e.hasNext();) {
                    String cellKey = (String) e.next();
                    IMSuite fullSuite = (IMSuite) trModel.get(cellKey);
                    int fullposSubs = fullSuite.getCompletePositiveSubjects(i);
                    int fullnegSubs = fullSuite.getCompleteNegativeSubjects(i);
                    double fullposScr = fullSuite.getCompletePositiveScore(i);
                    double fullnegScr = fullSuite.getCompleteNegativeScore(i);
                    if (tModel.containsKey(cellKey)) {
                        IMSuite testingSuite = (IMSuite) tModel.get(cellKey);
                        int pos_Subs;
                        int neg_Subs;
                        double pos_Scr;
                        double neg_Scr;
                        if (testingSuite.ExistIt(i, j)) {
                            pos_Subs = testingSuite.getPositiveSubjects(i, j);
                            neg_Subs = testingSuite.getNegativeSubjects(i, j);
                            pos_Scr = testingSuite.getPositiveScore(i, j);
                            neg_Scr = testingSuite.getNegativeScore(i, j);
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
                hc.SetStatistic(i, j, PublicData.TrainingAccuIdx, trAccu);
                hc.SetStatistic(i, j, PublicData.TestingAccuIdx, tAccu);
                if (currBestStatistic[i][j] < trAccu) {
                    currBestStatistic[i][j] = trAccu;
                    bestTestingAccu[j] = tAccu;
                    bestModel.put(new Integer(j), modelName);
                }
            }
        }
        hc.summarise();
        ModelHubs.put(modelName, hc);
    }

    private void ChrInt() {
        int c = 0;
        int idx = 0;
        hasMoreThanOneInterval = false;
        for (int i = 0; i < imp.ChromosomeNumber(); i++) {
            for (int j = 0; j < imp.IntervalNumberAtChromosome(i); j++) {
                if (c == SNPIndex[idx]) {
                    ChrInt[idx][0] = i; // chromosome
                    ChrInt[idx][1] = j; // interval
                    idx++;
                    if (idx > 0 && ChrInt[idx - 1][0] == i) {
                        hasMoreThanOneInterval = true;
                    }
                    if (idx == SNPIndex.length) {
                        break;
                    }
                }
                c++;    // code of combination of chromosome and interval
            }
            if (idx == SNPIndex.length) {
                break;
            }
        }
    }

    public void print() {
        List com = (List) comGenerator.get(order);
        for (Iterator e = com.iterator(); e.hasNext();) {
            String key = (String) e.next();
            if (!ModelHubs.containsKey(key)) {
                continue;
            }
            SNPIndex = ToolKit.StringToIntArray(key);
            ChrInt();
            for (int i = 0; i < ChrInt.length; i++) {
                for (int j = 0; j < ChrInt[i].length; j++) {
                    System.out.print(ChrInt[i][j] + ",");
                }
                System.out.print("\t");
            }
            System.out.println();
            HeteroCombination hc = ModelHubs.get(key);
            System.out.println(hc);
        }
    }

    public void printBestModel() {
        System.out.println("Best model");
        for (int i = 0; i < bestModel.size(); i++) {
            String key = bestModel.get(new Integer(i));
            SNPIndex = ToolKit.StringToIntArray(key);
            ChrInt();
            for (int j = 0; j < ChrInt.length; j++) {
                for (int k = 0; k < ChrInt[j].length; k++) {
                    System.out.print(ChrInt[j][k] + ",");
                }
                System.out.print("\t");
            }
            System.out.print(currBestStatistic[0][i] + " " + bestTestingAccu[i] + "\n");
        }
    }

    public void printWholeDataAccuracy() {
        System.out.println("whole data statistics");
        for (int i = 0; i < numTraits; i++) {
            HashMap<String, Double> t = (HashMap<String, Double>) wholedataAccuracy.get(i);
            List com = (List) comGenerator.get(order);
            for (Iterator e = com.iterator(); e.hasNext();) {
                String key = (String) e.next();
                if (!t.containsKey(key)) {
                    continue;
                }
                SNPIndex = ToolKit.StringToIntArray(key);
                ChrInt();
                for (int j = 0; j < ChrInt.length; j++) {
                    System.out.print(ChrInt[j][0] + "," + ChrInt[j][1] + "\t");
                }
                System.out.println((Double) t.get(key));
            }
        }
    }

    public ArrayList SignificantPoints(double threshold) {
        ArrayList sigPoints = new ArrayList();
        for (int i = 0; i < numTraits; i++) {
            HashMap<String, Double> t = (HashMap<String, Double>) wholedataAccuracy.get(i);
            List com = (List) comGenerator.get(order);
            for (Iterator e = com.iterator(); e.hasNext();) {
                String key = (String) e.next();
                StringBuffer loc = new StringBuffer();
                if (!t.containsKey(key)) {
                    continue;
                }
                if ((((Double) t.get(key)).doubleValue()) < threshold) {
                    continue;
                }
                SNPIndex = ToolKit.StringToIntArray(key);
                ChrInt();
//                for (int j = 0; j < ChrInt.length; j++) {
//                    System.out.print(ChrInt[j][0] + "," + ChrInt[j][1] + "\t");
//                }
//                System.out.println((Double) t.get(key));
                for (int j = 0; j < ChrInt.length; j++) {
                    loc.append(ChrInt[j][0] + PublicData.seperator + ChrInt[j][1] + PublicData.seperator);
                }
                sigPoints.add(loc.toString() + ((Double) t.get(key)) + PublicData.seperator + key);
            }
        }
        return sigPoints;
    }

    public double peakStatistic() {
        return peak;
    }

    public double getTopStatistic(int pheIdx, int idx) {
        double s = 0;
        List com = (List) comGenerator.get(order);
        for (Iterator e = com.iterator(); e.hasNext();) {
            String key = (String) e.next();
            SNPIndex = ToolKit.StringToIntArray(key);
            ChrInt();
            HeteroCombination hc = ModelHubs.get(key);
            s = s < hc.get(pheIdx, idx) ? hc.get(pheIdx, idx) : s;
        }
        return s;
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
