package mdr.moore.statistic;

public class MDRStatistic implements Comparable {
	private double testing_accuracy;
	private double training_accuracy;
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
	public void setTestingAccuracy(double ta) {
		testing_accuracy = ta;
	}
	
	public void setTrainingAccuracy(double tra) {
		training_accuracy = tra;
	}

	public double getTestingBalancedAccuracy() {
		return testing_balanced_accuracy;
	}
	
	public double getTrainingBalancedAccuracy() {
		return training_balanced_accuracy;
	}

	public double getTestingAccuracy() {
		return testing_accuracy;
	}
	
	public double getTrainingAccuracy() {
		return training_accuracy;
	}
	
	public int compareTo(Object mdrstat) {
		try {
		if(! (mdrstat instanceof MDRStatistic)) {
			throw new ClassCastException("Invalid MDRStatistic Object");
		}
		} catch (ClassCastException E) {
			E.printStackTrace(System.err);
		}
		
		double ta = ((MDRStatistic) mdrstat).getTestingAccuracy();
		if(testing_accuracy > ta) {
			return 1;
		} else if(testing_accuracy == ta){
			return 0;
		} else {
			return -1;
		}
	}
}
