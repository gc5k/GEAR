package family.mdr;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Set;
import java.util.Map.Entry;

import publicAccess.PublicData;

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

	private double[] statistic;

	private MDRStatistic mdr_stat;
	private String best_model;

	private int[] include_snp = null;
	private int[] exclude_snp = null;
	private int[] a;
	private int end = 0;
	private String includeComb = null;

	public HeteroLinearMergeSearch(DataFile dr, Subdivision sd, int[] inSNP, int[] exSNP, int num_marker) {
		super(dr, sd);
		int len = num_marker;
		if (inSNP != null) {
			include_snp = new int[inSNP.length];
			System.arraycopy(inSNP, 0, include_snp, 0, inSNP.length);
			len -= inSNP.length;
			StringBuffer s = new StringBuffer();
			for(int i = 0; i < inSNP.length; i++) {
				s.append(inSNP[i]);
				if( i != inSNP.length - 1) {
					s.append(PublicData.seperator);
				}
			}
			includeComb = s.toString();
		}
		if (exSNP != null) {
			exclude_snp = new int[exSNP.length];
			System.arraycopy(exSNP, 0, exclude_snp, 0, exSNP.length);
			len -= exSNP.length;
		}
		a = new int[len];
		int idx = 0;
		for (int i = 0; i < num_marker; i++) {
			boolean flag = true;
			if(include_snp != null) {
				for (int j = 0; j < include_snp.length; j++) {
					if (i == include_snp[j]) {
						flag = false;
						break;
					}
				}
			}
			if (!flag) continue;
			if(exclude_snp != null) {
				for (int j = 0; j < exclude_snp.length; j++) {
					if (i == exclude_snp[j]) {
						flag = false;
						break;
					}
				}
			}
			if (flag) {
				a[idx++] = i;
			}
		}
		end = a.length;
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

		int[] order = new int[or + 1];
		for (int i = 0; i <= or; i++) {
			order[i] = i - 1;
		}
		int count = 0;

		int k = or;
		
		int len = or;
		if (include_snp != null) {
			len -= include_snp.length;
		}
		boolean flag = true;
		StringBuilder combi;
		while (order[0] == -1) {
			if (flag) {
				combi = includeComb == null ? new StringBuilder() : new StringBuilder(includeComb);

				for (int i = 1; i <= len; i++) {
					combi.append(PublicData.seperator);
					combi.append(Integer.toString(a[order[i]]));
				}
				kernal(combi.toString());
				count++;
				flag = false;
			}
			order[k]++;
			if (order[k] == end) {
				order[k--] = 0;
				continue;
			}
			if (k < len) {
				order[k + 1] = order[k];
				k++;
				continue;
			}
			if (k == len) {
				flag = true;
			}
		}
	}

	private void kernal(String modelName) {
		ArrayList<DataFile.Subject> sample = data.getSample();

		cleanupTestingSet();

		model = new Combination();
		SNPIndex = ToolKit.StringToIntArray(modelName);
		linearSearch(sample);
		// mergeSearch(data.getSample(), null, 0);
		double[] t = calculateSingleBest(modelName);
		if (t[MDRConstant.TestingBalancedAccuIdx] > statistic[MDRConstant.TestingBalancedAccuIdx]) {
			statistic = t;
			best_model = modelName;
		}
		count++;
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
		if (mdr_stat == null) {
			mdr_stat = mdrstat;
			best_model = modelName;
		} else {
			if (mdr_stat.compareTo(mdrstat) < 0) {
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
		for (int i = 0; i < MDRConstant.NumStats; i++) {
			sb.append(MDRConstant.TestStatistic[i] + ", ");
		}
		sb.append(System.getProperty("line.separator"));
		for (Entry<String, MDRStatistic> entry : heteroresult.entrySet()) {
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
