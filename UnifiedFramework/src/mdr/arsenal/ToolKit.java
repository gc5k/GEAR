package mdr.arsenal;

import java.util.HashMap;
import java.util.Map.Entry;

import publicAccess.PublicData;

import util.NewIt;
import mdr.result.Cell;

/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class ToolKit {

	public static int[] StringToIntArray(String s) {
		String[] unit = s.split(PublicData.seperator);
		int[] value = new int[unit.length];
		for (int i = 0; i < unit.length; i++) {
			value[i] = Integer.parseInt(unit[i]);
		}
		return value;
	}

	public static int[] StringArrayTOIntArray(String[] s) {
		if (s == null) {
			return null;
		}
		int[] d = new int[s.length];
		for (int i = 0; i < d.length; i++) {
			d[i] = Integer.parseInt(s[i]);
		}
		return d;
	}

	public static String IntArrayToString(int[] array) {
		StringBuffer sb = new StringBuffer();
		for (int i = 0; i < array.length; i++) {
			Integer I = new Integer(array[i]);
			sb.append(I);
			if (i == (array.length - 1)) {
				continue;
			}
			sb.append(PublicData.seperator);
		}
		return sb.toString();
	}

	public static double getPMCC(double[] point1, double[] point2) {
		double suma = 0;
		double sumb = 0;
		double sumaSq = 0;
		double sumbSq = 0;
		double pSum = 0;
		int n = point1.length;
		for (int i = 0; i < point1.length; i++) {
			suma = suma + point1[i];
			sumb = sumb + point2[i];
			sumaSq = sumaSq + point1[i] * point1[i];
			sumbSq = sumbSq + point2[i] * point1[i];
			pSum = pSum + point1[i] * point2[i];
		}
		double numerator = pSum - suma * sumb / n;
		double denominator = Math.sqrt((sumaSq - suma * suma / n) * (sumbSq - sumb * sumb / n));
		return numerator / denominator;
	}

	public static double Mean(double[] vector) {
		double sum = 0;
		if (vector.length == 0) {
			return 0;
		}
		for (int i = 0; i < vector.length; i++) {
			sum += vector[i];
		}
		return sum / (double) vector.length;
	}

	public static double[] Mean(double[][] vector, boolean row) {
		double[] mean;
		double[][] v;
		if (!row) {
			v = transposeMatrix(vector);
		} else {
			v = vector;
		}
		mean = new double[v.length];
		for (int i = 0; i < mean.length; i++) {
			mean[i] += Mean(v[i]);
		}
		return mean;
	}

	public static double Variance(double[] vector) {
		double sum = 0, sumSquared = 0;
		if (vector.length <= 1) {
			return 0;
		}
		for (int i = 0; i < vector.length; i++) {
			sum += vector[i];
			sumSquared += (vector[i] * vector[i]);
		}
		double result = (sumSquared - (sum * sum / (double) vector.length)) / (double) (vector.length - 1);
		// We don't like negative variance
		if (result < 0) {
			return 0;
		} else {
			return result;
		}
	}

	public static double[] Variance(double[][] vector, boolean row) {
		double[] var;
		double[][] v;
		if (!row) {
			v = transposeMatrix(vector);
		} else {
			v = vector;
		}
		var = new double[v.length];
		for (int i = 0; i < var.length; i++) {
			var[i] = Variance(v[i]);
		}
		return var;
	}

	public static double[][] transposeMatrix(double[][] vector) {
		double[][] tm = new double[vector[0].length][vector.length];
		for (int i = 0; i < tm.length; i++) {
			for (int j = 0; j < tm[i].length; j++) {
				tm[i][j] = vector[j][i];
			}
		}
		return tm;
	}

	public static double BalancedAccuracy(HashMap<String, Cell> model) {
		try {
			if (model == null) {
				throw new ToolKitException("It is an empty model.");
			}
		} catch (ToolKitException E) {
			E.printStackTrace(System.err);
		}
		double accuracy = 0;
		double truePos = 0;
		double trueNeg = 0;
		double falsePos = 0;
		double falseNeg = 0;
		double unknown = 0;
		for (Entry<String, Cell> entry: model.entrySet()) {
			Cell cell = entry.getValue();
			if (cell.getStatus() == 1) {
				truePos += cell.getPositiveScore();
				falseNeg += Math.abs(cell.getNegativeScore());
			} else if (cell.getStatus() == 0) {
				trueNeg += Math.abs(cell.getNegativeScore());
				falsePos += cell.getPositiveScore();
			} else {
				unknown += Math.abs(cell.getNegativeScore()) + cell.getPositiveScore();
				// what's to do if for a unknown status;
			}
		}
		if( ((trueNeg + falsePos + unknown) == 0) || (truePos + falseNeg + unknown) == 0 ) {
			accuracy = 0;
		} else {
			accuracy = 0.5 * (truePos / (truePos + falseNeg + unknown/2) + trueNeg / (trueNeg + falsePos + unknown/2));
		}
		return accuracy;
	}

	public static double IMAccuracy(HashMap<String, Cell> model) throws ToolKitException {
		if (model == null) {
			throw new ToolKitException("It is an empty model.");
		}
		double accuracy = 0;
		double truePos = 0;
		double trueNeg = 0;
		double falsePos = 0;
		double falseNeg = 0;
		for (String key : model.keySet()) {
			Cell cell = (Cell) model.get(key);
			if (cell.getStatus() == 1) {
				truePos += Math.abs(cell.getPositiveScore() - cell.getExpectedPostiveScore())
						* Math.abs(cell.getPositiveScore() - cell.getExpectedPostiveScore());
				falseNeg += Math.abs((Math.abs(cell.getNegativeScore()) - Math.abs(cell.getExpectedNegativeScore())))
						* Math.abs((Math.abs(cell.getNegativeScore()) - Math.abs(cell.getExpectedNegativeScore())));
			} else if (cell.getStatus() == 0) {
				trueNeg += Math.abs(Math.abs(cell.getNegativeScore() - Math.abs(cell.getExpectedNegativeScore())))
						* Math.abs(Math.abs(cell.getNegativeScore() - Math.abs(cell.getExpectedNegativeScore())));
				falsePos += Math.abs(cell.getPositiveScore() - cell.getExpectedPostiveScore())
						* Math.abs(cell.getPositiveScore() - cell.getExpectedPostiveScore());
			} else {
				// what's going on if for a negative status;
			}
		}
		double denominator = truePos + trueNeg + falsePos + falseNeg;
		if (denominator == 0) {
			throw new ToolKitException("denominator is zero.");
		}
		accuracy = (truePos + trueNeg) / (denominator);
		return accuracy;
	}

	public static HashMap<String, Double> AccuracyLou(HashMap<String, Cell> model, HashMap<String, Double> QTLProbabilityMatrix)
			throws ToolKitException {
		HashMap<String, Double> posteriorProbability = NewIt.newHashMap();
		if (model == null) {
			throw new ToolKitException("It is an empty model.");
		}
		double weightH = 0;
		double weightL = 0;
		double truePos = 0;
		double trueNeg = 0;
		double falsePos = 0;
		double falseNeg = 0;
		for (String key : model.keySet()) {
			Cell cell = (Cell) model.get(key);
			if (cell.getStatus() == 1) {
				weightH += ((Double) QTLProbabilityMatrix.get(key)).doubleValue();
				truePos += cell.getPositiveScore();
				falseNeg += Math.abs(cell.getNegativeScore());
			} else if (cell.getStatus() == 0) {
				weightL += ((Double) QTLProbabilityMatrix.get(key)).doubleValue();
				trueNeg += Math.abs(cell.getNegativeScore());
				falsePos += cell.getPositiveScore();
			} else {
				// what's going on if for a negative status;
			}
		}
		for (String key : model.keySet()) {
			Cell cell = (Cell) model.get(key);
			if (cell.getStatus() == 1) {
				posteriorProbability.put(key, new Double(weightH));
			} else if (cell.getStatus() == 0) {
				posteriorProbability.put(key, new Double(weightL));
			} else {

			}
		}
		return posteriorProbability;
	}

	public static double CGBStatistic(HashMap<String, Cell> model) throws ToolKitException {
		double cgbStatistic = 0;

		return cgbStatistic;
	}

	public static int[] Sort(/* @non_null@ */double[] array) {
		int[] index = new int[array.length];
		array = (double[]) array.clone();
		for (int i = 0; i < index.length; i++) {
			index[i] = i;
			if (Double.isNaN(array[i])) {
				array[i] = Double.MAX_VALUE;
			}
		}
		QuickSort(array, index, 0, array.length - 1);
		return index;
	}

	private static void QuickSort(double[] array, int[] index, int left, int right) {
		if (left < right) {
			int middle = Partition(array, index, left, right);
			QuickSort(array, index, left, middle);
			QuickSort(array, index, middle + 1, right);
		}
	}

	private static int Partition(double[] array, int[] index, int l, int r) {
		double pivot = array[index[(l + r) / 2]];
		int help;
		while (l < r) {
			while ((array[index[l]] < pivot) && (l < r)) {
				l++;
			}
			while ((array[index[r]] > pivot) && (l < r)) {
				r--;
			}
			if (l < r) {
				help = index[l];
				index[l] = index[r];
				index[r] = help;
				l++;
				r--;
			}
		}
		if ((l == r) && (array[index[r]] > pivot)) {
			r--;
		}
		return r;
	}
}
