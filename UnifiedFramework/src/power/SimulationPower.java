package power;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map.Entry;

import family.mdr.result.MDRStatistic;

import util.NewIt;


/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public final class SimulationPower {
	private static String model = "0,3";
	private HashMap<String, MDRStatistic> result;
	private double[] threshold;
	private double p_001 = 0.01;
	private double t_001;
	private int p_001_idx = 0;
	private double p_005 = 0.05;
	private double t_005;
	private int p_005_idx = 0;
	public static int typeI_005;
	public static int typeI_001;

	public static int power_005;
	public static int power_001;

	private ArrayList<SimulationStatistic> result_005 = NewIt.newArrayList();
	private ArrayList<SimulationStatistic> result_001 = NewIt.newArrayList();

	public class SimulationStatistic {
		protected String model;
		private double TestingBalancedAccuracy;
		private double TrainingBalancedAccuracy;
		private double p_value;

		public SimulationStatistic(String m, double tra, double ta, double p) {
			model = m;
			TestingBalancedAccuracy = tra;
			TrainingBalancedAccuracy = ta;
			p_value = p;
		}

		public String toString() {
			NumberFormat number = new DecimalFormat("#.####");
			StringBuilder sb = new StringBuilder();
			sb.append(model + " " + number.format(TestingBalancedAccuracy) + " " + number.format(TrainingBalancedAccuracy) + " "
					+ number.format(p_value) + System.getProperty("line.separator"));
			return sb.toString();
		}
	}

	public SimulationPower(HashMap<String, MDRStatistic> r, double[] t) {
		result = r;
		threshold = t;
		Arrays.sort(threshold);
		p_005_idx = (int) (threshold.length * (1 - p_005));
		p_001_idx = (int) (threshold.length * (1 - p_001));
		t_005 = threshold[p_005_idx];
		t_001 = threshold[p_001_idx];
	}

	public void calculatePower() {

		for (Entry<String, MDRStatistic> entry : result.entrySet()) {
			String k = entry.getKey();
			MDRStatistic v = entry.getValue();
			double ta = v.getTestingBalancedAccuracy();
			if (ta > t_005) {
				result_005.add(new SimulationStatistic(k, v.getTestingBalancedAccuracy(), v.getTrainingBalancedAccuracy(), getP(p_005_idx, ta)));
			}
			if (ta > t_001) {
				result_001.add(new SimulationStatistic(k, v.getTestingBalancedAccuracy(), v.getTrainingBalancedAccuracy(), getP(p_001_idx, ta)));
			}
		}
		if (result_005.size() > 0) {
			typeI_005++;
			for (SimulationStatistic s : result_005) {
				if (s.model.compareTo(model) == 0) {
					power_005++;
				}
			}
		}
		if (result_001.size() > 0) {
			typeI_001++;
			for (SimulationStatistic s : result_001) {
				if (s.model.compareTo(model) == 0) {
					power_001++;
				}
			}
		}

	}

	private double getP(int idx, double t) {
		double p = idx;
		for (int i = idx, len = threshold.length; i < len; i++) {
			if (t > threshold[i])
				p++;
		}
		return 1 - p / threshold.length;
	}

	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(threshold.length + " rounds of permutation:");
		sb.append(System.getProperty("line.separator"));
		sb.append("significance level at alpha = 0.05, BTA threshold=" + String.format("%.4f", t_005));
		sb.append(System.getProperty("line.separator"));
		if (result_005.size() > 0) {
			for (SimulationStatistic sr : result_005) {
				sb.append(sr);
			}
		} else {
			sb.append("null");
			sb.append(System.getProperty("line.separator"));
		}
		sb.append("significance level at alpha = 0.01, BTA threshold=" + String.format("%.4f", t_001));
		sb.append(System.getProperty("line.separator"));
		if (result_001.size() > 0) {
			for (SimulationStatistic sr : result_001) {
				sb.append(sr);
			}
		} else {
			sb.append("null");
			sb.append(System.getProperty("line.separator"));
		}
		return sb.toString();
	}
}
