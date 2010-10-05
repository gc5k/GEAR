package mdr;

import java.util.ArrayList;

/**
 *
 * @author Guo-Bo Chen
 */
public class Cell {

    private int positiveSubjects;
    private int negativeSubjects;
    private int status;// 1 for high-risk; 0 for low-risk; -1 for empty
    private double positiveScore;
    private double negativeScore;
    private double expectedPosScore;
    private double expectedNegScore;

    public Cell(ArrayList subjects) {
    }

    public Cell() {
        status = -1;
    }

    public Cell(int posSubjects, int negSubjects, double posScore, double negScore, int s) {
        positiveSubjects = posSubjects;
        negativeSubjects = negSubjects;
        positiveScore = posScore;
        negativeScore = negScore;
        status = s;
    }

    public Cell(int posSubjects, int negSubjects, double posScore, double negScore, double ExpectedPosScore, double ExpectedNegScore, int s) {
        positiveSubjects = posSubjects;
        negativeSubjects = negSubjects;
        positiveScore = posScore;
        negativeScore = negScore;
        expectedPosScore = ExpectedPosScore;
        expectedNegScore = ExpectedNegScore;
        status = s;
    }
    
    public Cell(int affected, int unaffected) {
    }

    public double getPositiveSubjects() {
        return positiveSubjects;
    }

    public double getNegativeSubjects() {
        return negativeSubjects;
    }

    public double getPositiveScore() {
        return positiveScore;
    }

    public double getNegativeScore() {
        return negativeScore;
    }

    public double getExpectedPostiveScore() {
    	return expectedPosScore;
    }
    
    public double getExpectedNegativeScore() {
    	return expectedNegScore;
    }
    
    public int getStatus() {
        return status;
    }

    public String toString() {
        return "Cell status: " + status + ", Positive subjects: " + Integer.toString(positiveSubjects) + ", " + "Negative subjects: " + Integer.toString(negativeSubjects) + ", " + "Positive score: " + Double.toString(positiveScore) + ", " + "Negative score: " + Double.toString(negativeScore);
    }
}
