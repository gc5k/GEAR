package family.mdr;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.FDistributionImpl;
/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class MDRStatistic implements Comparable<MDRStatistic> {

	private double[] stats;
	private double Vt;
	private double Vx;
	private int N;
	private int Npos;
	private int Nneg;
	private double mean;
	private double FValue;
	private double PF;
	
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

	public void setVt(double vt) {
		Vt = vt;
	}
	
	public void setVx(double vx) {
		Vx = vx;
	}
	
	public void setN(int n) {
		N = n;
	}
	
	public void setMean(double m) {
		mean = m;
	}

	public void setNpos(int Np) {
		Npos = Np;
	}
	
	public void setNneg(int Nn) {
		Nneg = Nn;
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

	private void F() {
		double numerator = 0;
		double denominator = 0;
		int degreeNumerator = 0;
		int degreeDenominator = 0;

		if(Npos>0 && Nneg>0) {
			numerator = Vx / 2;
			denominator = Vt / (N-2);
			degreeNumerator = 2;
			degreeDenominator = N - 2;
			FValue = numerator / denominator;
		} else if(Npos > 0 && Nneg == 0) {
			numerator = Vx / 1;
			denominator = Vt / (N-1);
			degreeNumerator = 1;
			degreeDenominator = N - 1;
			FValue = numerator / denominator;
		} else if(Npos == 0 && Nneg > 0) {
			numerator = Vx / 1;
			denominator = Vt / (N-1);
			degreeNumerator = 1;
			degreeDenominator = N - 1;
			FValue = numerator / denominator;
		} else {
			degreeNumerator = 1;
			degreeDenominator = 1;
			FValue = 0;
		}
		FDistributionImpl F = new FDistributionImpl(degreeNumerator, degreeDenominator);

		try {
			PF = 1 - F.cumulativeProbability(FValue);
		} catch (MathException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		F();
		sb.append(String.format("%.4f", stats[0]) + ", " + String.format("%.4f", stats[1]) +  ", " + String.format("%.4f", Vt) +  ", " + String.format("%.4f", Vx) + ", " + N + ", " + String.format("%.4f", mean) + ", " + String.format("%.4f", FValue) + ", " + String.format("%.4E", PF));
		return sb.toString();
	}
}
