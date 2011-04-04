
package mdr;

import mdr.data.DataFile;
import java.util.ArrayList;
import java.util.Iterator;

import publicAccess.PublicData;
import util.NewIt;

/**
 *
 * @author Guo-Bo Chen
 */
public class Suite {
    protected ArrayList<DataFile.Subject> subjects;
    protected double[] threshold;
    protected int[] posSubjects;
    protected int[] negSubjects;
    protected double[] posScore;
    protected double[] negScore;
    protected int[] status;

    public Suite(Suite s) {
        subjects = s.getSubjects();
        initial(s.getNumTraits());
    }

    public Suite(ArrayList<DataFile.Subject> subs, int numTraits) {
        subjects = subs;
        initial(numTraits);
    }

    public Suite(int numTraits) {
        subjects = NewIt.newArrayList();
        initial(numTraits);
    }
    
    private void initial(int numTraits) {
        threshold = new double[numTraits];
        posSubjects = new int[numTraits];
        negSubjects = new int[numTraits];
        posScore = new double[numTraits];
        negScore = new double[numTraits];   
        status = new int[numTraits];
    }

    public void add(DataFile.Subject s) {
        subjects.add(s);
    }

    public int getNumTraits() {
        return posScore.length;
    }

    public Object get(int idx) {
        return subjects.get(idx);
    }

    public int getStatus(int idx) {
        return status[idx];
    }

    public ArrayList<DataFile.Subject> getSubjects() {
        return subjects;
    }

    public int size() {
        return subjects.size();
    }

    public void setStatus(int idx, int sts) {
        status[idx] = sts;
    }

    public void summarize(int idx, double offset) {
        posSubjects[idx] = 0;
        negSubjects[idx] = 0;
        posScore[idx] = 0;
        negScore[idx] = 0;
        for (Iterator e = subjects.iterator(); e.hasNext();) {
            DataFile.Subject s = (DataFile.Subject) e.next();
            Double sscore = s.getDoubleScore(idx);
            if (sscore != null) {
                double scr = sscore - offset;
                if ((scr) >= 0) {
                    posSubjects[idx]++;
                    posScore[idx] += scr;
                } else {
                    negSubjects[idx]++;
                    negScore[idx] += scr;
                }
            }
        }
    }

    public void ascertainment(int idx, double threshold) throws SuiteException {
        if (subjects == null) {
            throw new SuiteException("Empty ArrayList");
        }
        if (negScore[idx] == 0) {
            if (posScore[idx] == 0) {
                status[idx] = 0;
            } else {
                status[idx] = 1;
            }
        } else {
            if ((posScore[idx] / Math.abs(negScore[idx])) == threshold) {
                status[idx] = PublicData.tieValue;
            } else {
                status[idx] = (posScore[idx] / Math.abs(negScore[idx])) > threshold ? 1 : 0;
            }
        }
    }

    public double getNegativeScore(int idx) {
        return negScore[idx];
    }

    public double getPositiveScore(int idx) {
        return posScore[idx];
    }

    public int getNegativeSubjects(int idx) {
        return negSubjects[idx];
    }

    public int getPositiveSubjects(int idx) {
        return posSubjects[idx];
    }

    public String toString() {
        StringBuffer sb = new StringBuffer();
        for(DataFile.Subject sub : subjects ) {
            System.out.println(sub);
        }
        return sb.toString();
    }
    
    public void print() {
        for (int i = 0; i < posSubjects.length; i++) {
            System.out.println(posSubjects[i] + " " + negSubjects[i] + " " + posScore[i] + " " + negScore[i]);
        }
    }
}
