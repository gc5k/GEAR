
package mdr.result;

import mdr.data.DataFile;
import java.util.ArrayList;

import publicAccess.PublicData;
import util.NewIt;

/**
 *
 * @author Guo-Bo Chen
 */
public class Suite {
    protected ArrayList<DataFile.Subject> subjects;
    protected double threshold;
    protected int posSubjects;
    protected int negSubjects;
    protected double posScore;
    protected double negScore;
    protected int status;

    public Suite(Suite s) {
        subjects = s.getSubjects();
    }

    public Suite(ArrayList<DataFile.Subject> subs) {
        subjects = subs;
    }

    public Suite() {
        subjects = NewIt.newArrayList();
    }

    public void add(DataFile.Subject s) {
        subjects.add(s);
    }

    public Object get(int idx) {
        return subjects.get(idx);
    }

    public int getStatus() {
        return status;
    }

    public ArrayList<DataFile.Subject> getSubjects() {
        return subjects;
    }

    public int size() {
        return subjects.size();
    }

    public void setStatus(int sts) {
        status = sts;
    }

    public void summarize() {
        posSubjects = 0;
        negSubjects = 0;
        posScore = 0;
        negScore = 0;
        for (DataFile.Subject s: subjects) {

            double sscore = s.getSelectedScore();
            if (sscore != Double.NaN) {
                double scr = sscore;
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

    public void ascertainment() throws SuiteException {
        if (subjects == null) {
            throw new SuiteException("Empty ArrayList");
        }
        if (negScore == 0) {
            if (posScore == 0) {
                status = 0;
            } else {
                status = 1;
            }
        } else {
            if ((posScore / Math.abs(negScore)) == 1) {
                status = PublicData.tieValue;
            } else {
                status = (posScore / Math.abs(negScore)) > threshold ? 1 : 0;
            }
        }
    }

    public double getNegativeScore() {
        return negScore;
    }

    public double getPositiveScore() {
        return posScore;
    }

    public int getNegativeSubjects() {
        return negSubjects;
    }

    public int getPositiveSubjects() {
        return posSubjects;
    }

    public String toString() {
        StringBuffer sb = new StringBuffer();
        for(DataFile.Subject sub : subjects ) {
            System.out.println(sub);
        }
        return sb.toString();
    }
}
