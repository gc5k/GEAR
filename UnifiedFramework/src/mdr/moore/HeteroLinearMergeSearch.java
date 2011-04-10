package mdr.moore;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Set;


import mdr.MDRConstant;
import mdr.algorithm.Subdivision;
import mdr.arsenal.ToolKit;
import mdr.data.DataFile;

import mdr.moore.statistic.MDRStatistic;

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

	private HashMap<String, MDRStatistic> heteroresult = null;

	public HeteroLinearMergeSearch(DataFile dr, Subdivision sd, boolean ismooremdr) {
		super(dr, sd, ismooremdr);
	}

	public void search(int or, ArrayList<String> modelspace) {
		order = or;
		bestKFold = new BestKFoldCVResult(or, subdivision.getInterval());
		count = 0;
		Arrays.fill(currBestStats, 0);
		heteroresult = NewIt.newHashMap();

		for (String modelName : modelspace) {
			for (Combination testingModel : cvTestingSet) {
				testingModel.clear();
			}
			model = new Combination();
			SNPIndex = ToolKit.StringToIntArray(modelName);
			mergeSearch(data.getSample(), null, 0);
			calculateSingleBest(modelName);
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
		return mean;
	}

	public String toString() {
		StringBuilder sb = new StringBuilder();
		for(String m : heteroresult.keySet()) {
			sb.append(m);
			sb.append(" ");
			sb.append(heteroresult.get(m));
			sb.append("\n");
		}
		return sb.toString();
	}
}
