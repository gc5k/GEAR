package mdr.moore.statistic;

public class MDRStatistic implements Comparable<MDRStatistic> {
	private double testing_balanced_accuracy;
	private double training_balanced_accuracy;
	
	public MDRStatistic() {
		
	}
	
	public void setTestingBalancedAccuracy(double tba) {
		testing_balanced_accuracy = tba;
	}
	
	public void setTrainingBalancedAccuracy(double trba) {
		training_balanced_accuracy = trba;
	}

	public double getTestingBalancedAccuracy() {
		return testing_balanced_accuracy;
	}

	public double getTrainingBalancedAccuracy() {
		return training_balanced_accuracy;
	}
	
	public int compareTo(MDRStatistic mdrstat) {
		double ta = mdrstat.getTestingBalancedAccuracy();
		if(testing_balanced_accuracy > ta) {
			return 1;
		} else if(testing_balanced_accuracy == ta){
			return 0;
		} else {
			return -1;
		}
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(testing_balanced_accuracy + " " + training_balanced_accuracy);
		return sb.toString();
	}
}
