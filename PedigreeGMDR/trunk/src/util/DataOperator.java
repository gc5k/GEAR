package util;

import javastat.util.BasicStatistics;

public class DataOperator {
	
	public static double[][] Normalize(double[][] d, boolean center, boolean scale) {
		double[][] data = new double[d.length][];
		BasicStatistics bs = new BasicStatistics();
		for(int i = 0; i < d.length; i++) {
			double mean = center ? bs.mean(d[i]):0.0;
			double sd = scale ? Math.sqrt(bs.variance(d[i])):1.0;
			data[i] = new double[d[i].length];
			for(int j = 0; j < d[i].length; j++) {
				data[i][j] = (d[i][j] - mean)/sd;
			}
		}
		return data;
	}
	
	public static double[] Normalize(double[] d, boolean center, boolean scale) {
		double[] data = new double[d.length];
		BasicStatistics bs = new BasicStatistics();
		for(int i = 0; i < d.length; i++) {
			double mean = center ? bs.mean(d):0.0;
			double sd = scale ? Math.sqrt(bs.variance(d)):1.0;
			for(int j = 0; j < d.length; j++) {
				data[i] = (d[i] - mean)/sd;
			}
		}
		return data;
	}
}
