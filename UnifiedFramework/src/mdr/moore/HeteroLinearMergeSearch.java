package mdr.moore;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Set;
import java.util.Map.Entry;

import family.mdr.MDRStatistic;

import mdr.MDRConstant;
import mdr.algorithm.Subdivision;
import mdr.arsenal.ToolKit;
import mdr.data.DataFile;


import mdr.result.BestKFoldCVResult;
import mdr.result.Cell;
import mdr.result.Combination;
import mdr.result.OneCVSet;
import mdr.result.Suite;

import util.NewIt;

/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */

public class HeteroLinearMergeSearch extends LinearMergeSearch {

	private double[] statistic;

	private MDRStatistic mdr_stat;
	private String best_model;

	public HeteroLinearMergeSearch(DataFile dr, Subdivision sd) {
		super(dr, sd);
	}

	public void search(int or, ArrayList<String> modelspace) {
		statistic = new double[MDRConstant.NumStats];
		best_model = null;
		mdr_stat = null;

		order = or;
		bestKFold = new BestKFoldCVResult(or, subdivision.getInterval());
		count = 0;
		Arrays.fill(currBestStats, 0);

		heteroresult = NewIt.newHashMap();

		ArrayList<DataFile.Subject> sample = data.getSample();
		for (String modelName : modelspace) {
			for (Combination testingModel : cvTestingSet) {
				testingModel.clear();
			}
			model = new Combination();
			SNPIndex = ToolKit.StringToIntArray(modelName);
			linearSearch(sample);
//			mergeSearch(data.getSample(), null, 0);
			double[] t = calculateSingleBest(modelName);
			if(t[MDRConstant.TestingBalancedAccuIdx] > statistic[MDRConstant.TestingBalancedAccuIdx]) {
				statistic = t;
				best_model = modelName;
			}
			count++;
		}
	}

	public double[] calculateSingleBest(String modelName) {

		double[] mean = new double[MDRConstant.NumStats];
		for (int j = 0; j < subdivision.getInterval(); j++) {
			OneCVSet cvSet = new OneCVSet(j, modelName);
			Combination testingModels = (Combination) cvTestingSet.get(j);

			Cell trCell;
			Cell tCell;
			Set<String> cellKeys = model.keySet();
			double[] trStatus = new double[cellKeys.size()];
			double[] tStatus = new double[cellKeys.size()];
			int idx = 0;
			for (String cellKey : cellKeys) {
				int tr_status = -1;
				int t_status = -1;
				Suite fullSuite = (Suite) model.get(cellKey);
				int fullposSubs = fullSuite.getPositiveSubjects();
				int fullnegSubs = fullSuite.getNegativeSubjects();
				double fullposScr = fullSuite.getPositiveScore();
				double fullnegScr = fullSuite.getNegativeScore();
				if (testingModels.containsKey(cellKey)) {
					Suite testingSuite = (Suite) testingModels.get(cellKey);
					int pos_Subs = testingSuite.getPositiveSubjects();
					int neg_Subs = testingSuite.getNegativeSubjects();
					double pos_Scr = testingSuite.getPositiveScore();
					double neg_Scr = testingSuite.getNegativeScore();
					tr_status = ToolKit.Acertainment(fullposScr - pos_Scr, fullnegScr - neg_Scr);
					t_status = ToolKit.Acertainment(pos_Scr, neg_Scr);
					trCell = new Cell(fullposSubs - pos_Subs, fullnegSubs - neg_Subs, fullposScr - pos_Scr, fullnegScr - neg_Scr, tr_status);
					tCell = new Cell(pos_Subs, neg_Subs, pos_Scr, neg_Scr, tr_status);
				} else {
					tr_status = ToolKit.Acertainment(fullposScr, fullnegScr);
					trCell = new Cell(fullposSubs, fullnegSubs, fullposScr, fullnegScr, tr_status);
					tCell = new Cell(0, 0, 0, 0, -1);
				}
				cvSet.addTrainingModel(cellKey, trCell);
				cvSet.addTestingModel(cellKey, tCell);
				trStatus[idx] = tr_status;
				tStatus[idx] = t_status;
				idx++;
			}
			double trAccu = 0;
			double tAccu = 0;
			trAccu = ToolKit.BalancedAccuracy(cvSet.getTrainingSubdivision());
			mean[MDRConstant.TrainingBalancedAccuIdx] += trAccu;

			tAccu = ToolKit.BalancedAccuracy(cvSet.getTestingSubdivision());
			mean[MDRConstant.TestingBalancedAccuIdx] += tAccu;
			cvSet.setStatistic(MDRConstant.TrainingBalancedAccuIdx, trAccu);
			cvSet.setStatistic(MDRConstant.TestingBalancedAccuIdx, tAccu);
		}
		mean[MDRConstant.TrainingBalancedAccuIdx] /= subdivision.getInterval();
		mean[MDRConstant.TestingBalancedAccuIdx] /= subdivision.getInterval();
		MDRStatistic mdrstat = new MDRStatistic();
		heteroresult.put(modelName, mdrstat);
		mdrstat.setTrainingBalancedAccuracy(mean[MDRConstant.TrainingBalancedAccuIdx]);
		mdrstat.setTestingBalancedAccuracy(mean[MDRConstant.TestingBalancedAccuIdx]);
		if(mdr_stat == null) {
			mdr_stat = mdrstat;
			best_model = modelName;
		} else {
			if(mdr_stat.compareTo(mdrstat) < 0) {
				best_model = modelName;
				mdr_stat = mdrstat;
			}
		}
		return mean;
	}

	public String getBestModelKey() {
		return best_model;
	}

	public double[] getStats() {
		return statistic;
	}

	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("model, ");
		for(int i = 0; i < MDRConstant.NumStats; i++) {
			sb.append(MDRConstant.TestStatistic[i] + ", ");
		}
		sb.append(System.getProperty("line.separator"));
		for(Entry<String, MDRStatistic> entry : heteroresult.entrySet()) {
			sb.append(entry.getKey());
			sb.append(" ");
			sb.append(entry.getValue());
			sb.append(System.getProperty("line.separator"));
		}
		sb.append("===================");
		sb.append(System.getProperty("line.separator"));
		sb.append("best model: " + best_model);
		sb.append(System.getProperty("line.separator"));
		sb.append(mdr_stat);
		return sb.toString();
	}
}
