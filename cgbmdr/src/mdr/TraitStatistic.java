package mdr;

import mdr.data.DataFile;
import java.util.ArrayList;
import java.util.Iterator;
import im.population.IMPopulation;

/**
 *
 * @author Guo-Bo Chen
 */
public class TraitStatistic {

    private double posScore = 0;
    private int posSubjects = 0;
    private double negScore = 0;
    private int negSubjects = 0;
    private double offset = 0;
    private int scrIdx;

    public TraitStatistic(int sI, ArrayList sam) {
        scrIdx = sI;
        double sum = 0;
        int count = 0;
        for (Iterator e = sam.iterator(); e.hasNext();) {
            DataFile.Subject s = (DataFile.Subject) e.next();
            Double sscore = s.getDoubleScore(scrIdx);
            if (sscore != null) {
                sum += sscore;
                count++;
            }
        }

        offset = sum / count;
        offsetting(offset, sI, sam);
    }

    public TraitStatistic(double os, int sI, ArrayList sam) {
        offset = os;
        scrIdx = sI;
        offsetting(os, sI, sam);
    }

    public TraitStatistic() {
    }

    /**
     * helper method, offsetting the scores of instances
     *
     * @param os    offset value
     * @param sI    score index
     * @param sam   list of instances
     */
    private void offsetting(double os, int sI, ArrayList sam) {
        for (Iterator e = sam.iterator(); e.hasNext();) {
            DataFile.Subject s = (DataFile.Subject) e.next();
            Double sscore = s.getDoubleScore(sI);
            if (sscore != null) {
                double scr = sscore - os;
                if ((scr) >= 0) {
                    posSubjects++;
                    posScore += scr;
                } else {
                    negSubjects++;
                    negScore += scr;
                }
            }
        }
    }

    public void DefaultScore(IMPopulation imp, int sI) {
        double mean = CalculateMU(imp, sI);
        for (int i = 0; i < imp.IndividualNumber(); i++) {
            if(!imp.PhenotypeExist(i, sI)) {
                continue;
            }
            double scr = imp.PhenotypeAt(i, sI) - mean;
            if (scr > 0) {
                posSubjects++;
                posScore += scr;
            } else {
                negSubjects++;
                negScore += scr;
            }
        }
    }

    public double CalculateMU(IMPopulation imp, int sI) {
        double sum = 0;
        int c = 0;
        for (int i = 0; i < imp.IndividualNumber(); i++) {
            if(!imp.PhenotypeExist(i, sI)) {
                c++;
            } else {
                sum += imp.PhenotypeAt(i, sI);
            }
        }
        return sum / (imp.IndividualNumber() - c);
    }

    public double getTraitPositiveScore() {
        return posScore;
    }

    public double getTraitNegativeScore() {
        return negScore;
    }

    public int getTraitPositiveSubjects() {
        return posSubjects;
    }

    public int getTraitNegativeSubjects() {
        return negSubjects;
    }

    public double getOffset() {
        return offset;
    }
}
