package mdr.result;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import mdr.MDRConstant;

import publicAccess.PublicData;

import util.NewIt;

/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class BestKFoldCVResult {

	private static final long serialVersionUID = 1L;
	int order;
    private ArrayList<OneCVSet> cvset = NewIt.newArrayList();
	private double[] KFoldMean; 
	private double[] bestModelStatistics;
	private HashMap<String, Integer> modelCount = NewIt.newHashMap();
	private String bestModel;
	private int cvConsistency;

	public BestKFoldCVResult(int or, int interval) {
		order = or;
		cvset.ensureCapacity(interval);
	}

	public String getModelAtCV(int i) {
		OneCVSet cvSet = get(i);
		return cvSet.getModel();
	}

	public String getBestModel() {
		return bestModel;
	}

	public int getCVC() {
		return cvConsistency;
	}

	public void add(OneCVSet cv) {
		cvset.add(cv);
	}

	public OneCVSet get(int i) {
		return cvset.get(i);
	}

	public OneCVSet set(int j, OneCVSet cvSet) {
		return cvset.set(j, cvSet);
	}

	public int getOrder() {
		return order;
	}

	public double[] getKFoldStats() {
		return KFoldMean;
	}

	public void summarise() {
		KFoldStatistic();
		
		bestModelStatistics = new double[MDRConstant.NumStats];

		Integer v = Collections.max(modelCount.values());
		cvConsistency = v.intValue();
		HashMap<String, Integer> majorModel = NewIt.newHashMap();
		ArrayList<String> majorKey = NewIt.newArrayList();

		for (String key : modelCount.keySet()) {
			if (v == modelCount.get(key)) {
				majorModel.put(key, v);
				majorKey.add(key);
			}
		}

		double[][] _bestModelStats = new double[majorKey.size()][MDRConstant.NumStats];
		for (int i = 0; i < cvset.size(); i++) {
			OneCVSet cvSet = get(i);
			String key = cvSet.getModel();
			if (!majorKey.contains(key)) {
				continue;
			}
			int idx = majorKey.indexOf(key);
			for (int j = 0; j < MDRConstant.NumStats; j++) {
				_bestModelStats[idx][j] += cvSet.getStatistic(j);
			}
		}

		double[] bigStats = new double[MDRConstant.NumStats];
		int idx = 0;
		for (String key : majorKey) {
			for (int j = 0; j < _bestModelStats[idx].length; j++) {
				_bestModelStats[idx][j] /= majorModel.get(key).intValue();
				if (bigStats[j] < _bestModelStats[idx][j]) bigStats[j] = _bestModelStats[idx][j];
			}
			idx++;
		}

		for (int i = 0; i < MDRConstant.NumStats; i++) {
			int c = 0;
			for (int j = 0; j < _bestModelStats.length; j++) {
				if ((bigStats[i] - _bestModelStats[j][i]) < PublicData.epsilon) {
					c++;
					bestModel = majorKey.get(j);
				}
			}
			if (c == 1) {
				break;
			}
		}
		int bestidx = majorKey.indexOf(bestModel);

		System.arraycopy(_bestModelStats[bestidx], 0, bestModelStatistics, 0, _bestModelStats[bestidx].length);
	}

	/**
	 * Calculate the test statistic based on k fold results
	 * The statistics may be biased if there is heterogeneity 
	 */
	private void KFoldStatistic() {
		KFoldMean = new double[MDRConstant.NumStats];

		for (int i = 0; i < cvset.size(); i++) {
			OneCVSet cvSet = (OneCVSet) get(i);
			for (int j = 0; j < MDRConstant.NumStats; j++) {
				KFoldMean[j] += cvSet.getStatistic(j);
			}

			String key = cvSet.getModel();
			Integer count = modelCount.get(key);
			if (count == null) {
				modelCount.put(key, new Integer(1));
			} else {
				int v = count.intValue();
				modelCount.put(key, new Integer(++v));
			}
		}

		for (int i = 0; i < MDRConstant.NumStats; i++) {
			KFoldMean[i] /= cvset.size();
		}
	}
}
