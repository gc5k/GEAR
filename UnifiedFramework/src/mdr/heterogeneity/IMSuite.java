package mdr.heterogeneity;

import java.util.ArrayList;
import java.util.HashMap;
/**
 *
 * @author Guo-Bo Chen
 */
public class IMSuite {

    double[][] posScore;
    double[][] negScore;
    double[][] posProbability;
    double[][] negProbability;
    int[][] posSubjects;
    int[][] negSubjects;
    ArrayList<HashMap> IndexMap;
    double[] posScr;
    double[] negScr;
    double[] posProb;
    double[] negProb;
    int[] posSubs;
    int[] negSubs;

    public IMSuite(int numTraits, int interval) {
        posScore = new double[numTraits][interval];
        negScore = new double[numTraits][interval];
        posProbability = new double[numTraits][interval];
        negProbability = new double[numTraits][interval];
        posSubjects = new int[numTraits][interval];
        negSubjects = new int[numTraits][interval];
        IndexMap = new ArrayList();
    }

    public IMSuite(int numTraits) {
        posScr = new double[numTraits];
        negScr = new double[numTraits];
        posProb = new double[numTraits];
        negProb = new double[numTraits];
        posSubs = new int[numTraits];
        negSubs = new int[numTraits];
    }

    public void addScore(int traitIdx, int interval, double score, double probability) {
        if (score > 0) {
            posScore[traitIdx][interval] += score;
            posProbability[traitIdx][interval] += probability;
            posSubjects[traitIdx][interval]++;
        } else {
            negScore[traitIdx][interval] += score;
            negProbability[traitIdx][interval] += probability;
            negSubjects[traitIdx][interval]++;
        }
    }

    public void addScore(int traitIdx, double score, double probability) {
        if (score > 0) {
            posScr[traitIdx] += score;
            posProb[traitIdx] += probability;
            posSubs[traitIdx]++;
        } else {
            negScr[traitIdx] += score;
            negProb[traitIdx] += probability;
            negSubs[traitIdx]++;
        }
    }

    public double getCompleteNegativeScore(int idx) {
        return negScr[idx];
    }

    public double getCompletePositiveScore(int idx) {
        return posScr[idx];
    }

    public int getCompleteNegativeSubjects(int idx) {
        return negSubs[idx];
    }

    public int getCompletePositiveSubjects(int idx) {
        return posSubs[idx];
    }

    public double getNegativeScore(int traitIdx, int intvl) {
        return negScore[traitIdx][intvl];
    }

    public double getPositiveScore(int traitIdx, int intvl) {
        return posScore[traitIdx][intvl];
    }

    public int getNegativeSubjects(int traitIdx, int intvl) {
        return negSubjects[traitIdx][intvl];
    }

    public int getPositiveSubjects(int traitIdx, int intvl) {
        return posSubjects[traitIdx][intvl];
    }

    public double getPositiveProbability(int traitIdx) {
    	return posProb[traitIdx];
    }

    public double getNegtiveProbability(int traitIdx) {
    	return negProb[traitIdx];
    }
    
    public boolean ExistIt(int traitIdx, int intvl) {
        boolean flag = true;
        if(posSubjects[traitIdx][intvl] == 0 && negSubjects[traitIdx][intvl] == 0) {
            flag = false;
        }
        return flag;
    }
}
