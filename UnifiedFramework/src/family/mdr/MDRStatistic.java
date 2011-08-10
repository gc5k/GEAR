package family.mdr;
/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class MDRStatistic implements Comparable<MDRStatistic> {

	private double[] stats;  
	//stats[0] for testing balanced accuracy; stats[1] for training balanced accuracy
	public MDRStatistic() {
		stats = new double[2];
	}
	
	public MDRStatistic(double tba, double trba) {
		stats = new double[2];
		stats[0] = tba;
		stats[1] = trba;
	}
	
	public void setTestingBalancedAccuracy(double tba) {
		stats[0] = tba;
	}
	
	public void setTrainingBalancedAccuracy(double trba) {
		stats[1] = trba;
	}

	public double getTestingBalancedAccuracy() {
		return stats[0];
	}

	public double getTrainingBalancedAccuracy() {
		return stats[1];
	}
	
	public double[] getStats() {
		return stats;
	}

	public int compareTo(MDRStatistic mdrstat) {
		double ta = mdrstat.getTestingBalancedAccuracy();
		if(stats[0] > ta) {
			return 1;
		} else if(stats[0] == ta){
			return 0;
		} else {
			return -1;
		}
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(String.format("%.4f", stats[0]) + ", " + String.format("%.4f", stats[1]));
		return sb.toString();
	}
}
