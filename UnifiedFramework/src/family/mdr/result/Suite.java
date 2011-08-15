
package family.mdr.result;

import family.mdr.MDRConstant;

import java.util.ArrayList;

import family.mdr.data.PersonIndex;

import util.NewIt;

/**
 *
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class Suite {
    protected ArrayList<PersonIndex> subjects;
    protected double threshold = 1;
    protected int posSubjects;
    protected int negSubjects;
    protected double posScore;
    protected double negScore;
    protected int status;

    public Suite(Suite s) {
        subjects = s.getSubjects();
    }

    public Suite(ArrayList<PersonIndex> subs) {
        subjects = subs;
    }

    public Suite() {
        subjects = NewIt.newArrayList();
    }

    public void add(PersonIndex s) {
        subjects.add(s);
    }

    public Object get(int idx) {
        return subjects.get(idx);
    }

    public int getStatus() {
        return status;
    }

    public ArrayList<PersonIndex> getSubjects() {
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
        for (PersonIndex s: subjects) {

            double sscore = s.getScore();
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

    public void ascertainment() {
    	status = Ascertainment(posScore, negScore);
    }

    public static int Ascertainment(double ps, double ns) {
		int s = 0;
		if (ns == 0) {
			if (ps == 0) {
				s = -1;
			} else {
				s = 1;
			}
		} else {
			if ((ps / Math.abs(ns)) == 1) {
				s = MDRConstant.tieValue;
			} else {
				s = (ps / Math.abs(ns)) > 1 ? 1 : 0;
			}
		}
		return s;
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
    	ascertainment();
        StringBuffer sb = new StringBuffer();
        if(status == 1) {
        	sb.append("H, ");
        } else if( status == 0) {
        	sb.append("L, ");
        } else {
        	sb.append("NA, ");
        }
        sb.append(String.format("%.2f", posScore) + "," + String.format("%d", posSubjects) + "," + String.format("%.2f", negScore) + "," + String.format("%d", -1*negSubjects));
//        for(DataFile.Subject sub : subjects ) {
//            System.out.println(sub);
//        }
        return sb.toString();
    }
}
