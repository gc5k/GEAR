package family.report;
import java.util.ArrayList;
import java.util.HashMap;

import publicAccess.PublicData;
import publicAccess.ToolKit;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;

public class Report {
	
	protected StringBuffer sum;
	protected ArrayList<ArrayList<Double>> test_statistics;
	protected ArrayList<ArrayList<Double>> null_test_statistics;

	protected String curr_model;
	public Report() {
		sum = new StringBuffer();
		sum.append("model\tTA\tZ\tP-value\tTrA\tZ-score\tP-value" + System.getProperty("line.separator"));
		test_statistics = new ArrayList();
	}

	public void NewRound(String m, boolean create_null) {
		curr_model = m;
		if(create_null) {
			null_test_statistics = new ArrayList();
		}
	}

	public void Add_test_statistic(double[] stat) {
		ArrayList<Double> ts = new ArrayList();
		for(int i = 0; i < stat.length; i++) {
			ts.add(new Double(stat[i]));
		}
		test_statistics.add(ts);
	}

	public void Add_null_test_statistic(double[] stat) {
		ArrayList<Double> ts = new ArrayList();
		for(int i = 0; i < stat.length; i++) {
			ts.add(new Double(stat[i]));
		}
		null_test_statistics.add(ts);
	}

	public void RoundSummary() {
		double[][] stats = new double[null_test_statistics.size()][PublicData.NumOfStatistics];
		for (int i = 0; i < stats.length; i++) {
			ArrayList t = null_test_statistics.get(i);
			for (int j = 0; j < stats[i].length; j++) {
				stats[i][j] = ((Double) t.get(j)).doubleValue();
			}
		}
		double[] mean = ToolKit.Mean(stats, false);
		double[] var = ToolKit.Variance(stats, false);
		ArrayList<Double> curr_stats = test_statistics.get(test_statistics.size() - 1);
		NormalDistributionImpl ND = new NormalDistributionImpl();
		double p = 0;
		sum.append(curr_model + "\t");
		for(int i = 0; i < curr_stats.size(); i++) {
			sum.append(curr_stats.get(i).toString() + "\t");
			double z = (((Double)curr_stats.get(i)).doubleValue() - mean[i])/Math.sqrt(var[i]);
			sum.append(z + "\t");
			try {
				p = 1 - ND.cumulativeProbability(z);
			} catch (MathException E) {
				E.printStackTrace(System.err);
			}
			sum.append(p + "\t");
		}
		sum.append(System.getProperty("line.separator"));
	}
	
	public String toString() {
		return sum.toString();
	}
}
