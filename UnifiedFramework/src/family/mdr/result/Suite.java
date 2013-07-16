
package family.mdr.result;


import java.text.DecimalFormat;
import java.util.ArrayList;

import admixture.parameter.Parameter;
import family.mdr.data.PersonIndex;

import util.NewIt;

/**
 *
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class Suite {
    protected ArrayList<PersonIndex> subjects;
    protected static double threshold = 1;
    protected int posSubjects;
    protected int negSubjects;
    protected double posScore;
    protected double negScore;
    protected double meanScore;
    protected int status = -1;

    public Suite(Suite s) {
        subjects = s.getSubjects();
    }

    public Suite(ArrayList<PersonIndex> subs) {
        subjects = subs;
    }

    public Suite() {
        subjects = NewIt.newArrayList();
    }

    public static void setThreshold(double T) {
    	threshold = T;
    }

    public static double getThreshold() {
    	return threshold;
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
                if ((sscore) >= 0) {
                    posSubjects++;
                    posScore += sscore;
                } else {
                    negSubjects++;
                    negScore += sscore;
                }
            }

        }
        meanScore = (posScore + negScore)/(posSubjects + negSubjects);
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
			if (ps == threshold * Math.abs(ns)) {
				s = Parameter.tie;
			} else {
				s = ps > threshold  * Math.abs(ns)? 1 : 0;
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
    
    public double getMeanScore() {
    	return meanScore;
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
        
		DecimalFormat fmt = new DecimalFormat("#.###E0");
		
//        sb.append(String.format("%.2f", posScore) + "," + String.format("%d", posSubjects) + "," + String.format("%.2f", negScore) + "," + String.format("%d", -1*negSubjects) + "," + String.format("%.4f", meanScore));
        sb.append(fmt.format(posScore) + "," + fmt.format(posSubjects) + "," + fmt.format(negScore) + "," + fmt.format(negSubjects));        
        return sb.toString();
    }
}
