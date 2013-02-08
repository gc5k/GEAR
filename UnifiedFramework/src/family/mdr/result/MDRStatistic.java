package family.mdr.result;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.FDistributionImpl;

import admixture.parameter.Parameter;
/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class MDRStatistic implements Comparable<MDRStatistic> {

	private double[] stats;
	private double meanT = Double.NaN;
	private double seT = Double.NaN;
	private double pT = Double.NaN;
	private double Vt;
	private double Vx;
	private int N;
	private int Npos;
	private int Nneg;
	private double mean;
	private double FValue;
	private int degreeDenominator;
	private int degreeNumerator;
	private double PF;
	
	private double TrainingP;
	private double TestingP;
	private double TruncatedFisherOneTailP;
	private double TruncatedFisherTwoTailP;
	//stats[0] for testing balanced accuracy; stats[1] for training balanced accuracy
	public MDRStatistic() {
		stats = new double[2];
	}
	
	public MDRStatistic(double tba, double trba) {
		stats = new double[2];
		stats[0] = tba;
		stats[1] = trba;
	}
	
	public void setTrainingPValue(double p) {
		TrainingP = p;
	}
	
	public double getTrainingPValue() {
		return TrainingP;
	}
	
	public void setTestingPValue(double p) {
		TestingP = p;
	}
	
	public double getTestingPValue() {
		return TestingP;
	}

	public void setTestingBalancedAccuracy(double tba) {
		stats[0] = tba;
	}
	
	public void setTestingBalancedAccuracyMeanT(double meanT) {
		this.meanT = meanT;
	}
	
	public void setTestingBalancedAccuracyseT(double seT) {
		this.seT = seT;
	}

	public void setTestingBalancedAccuracyPT(double pT) {
		this.pT = pT;
	}
	
	public double getTestingBalancedPT() {
		return this.pT;
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
	
	public double getVc() {
		return Vx/Vt;
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

	public void setTruncatedFisherOneTailP(double p) {
		TruncatedFisherOneTailP = p;
	}

	public void setTruncatedFisherTwoTailP(double p) {
		TruncatedFisherTwoTailP = p;
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

		if(Npos>0 && Nneg>0) {
			numerator = Vx / 2;
			denominator = Vt / (N-2);
			degreeNumerator = 1;
			degreeDenominator = N - 2;
			FValue = numerator / denominator;
//		} else if (Npos > 0 && Nneg == 0) {
//			numerator = Vx / 1;
//			denominator = Vt / (N - 1);
//			degreeNumerator = 1;
//			degreeDenominator = N - 1;
//			FValue = numerator / denominator;
//		} else if (Npos == 0 && Nneg > 0) {
//			numerator = Vx / 1;
//			denominator = Vt / (N - 1);
//			degreeNumerator = 1;
//			degreeDenominator = N - 1;
//			FValue = numerator / denominator;
		} else {
			degreeNumerator = 1;
			degreeDenominator = 1;
			FValue = 0;
		}
		FDistributionImpl F = new FDistributionImpl(degreeNumerator, degreeDenominator);

		try {
			PF = 1 - F.cumulativeProbability(FValue);
		} catch (MathException e) {
			e.printStackTrace();
		}
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		F();
		sb.append(N + ", ");
		sb.append(String.format("%.4f", Vx/Vt) + ", ");
		sb.append(String.format("%.4f", Vx));
		sb.append(", ");
		sb.append(String.format("%.4f", Vt));
		sb.append(", ");
//		sb.append(String.format("%.4f", FValue));
//		sb.append(", ");
//		sb.append(String.format("%.4E", PF));
//		sb.append(", "+degreeNumerator +", "+degreeDenominator+", ");
		sb.append(String.format("%.4f", stats[1]) + ", ");
		sb.append(String.format("%.4f", stats[0]) + ", ");
		if (Parameter.permFlag) {
			sb.append(String.format("%.4f", meanT) + ", " 
				+ String.format("%.4f", seT)+ ", " + String.format("%.4E", pT) + ", ");
		}
		

		return sb.toString();
	}
}
