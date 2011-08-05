package family.mdr;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Set;
import java.util.Map.Entry;

import admixture.parameter.Parameter;

import mdr.MDRConstant;
import mdr.algorithm.Subdivision;
import mdr.arsenal.ToolKit;
import mdr.data.DataFile;

import mdr.result.BestKFoldCVResult;
import mdr.result.Cell;
import mdr.result.Combination;
import mdr.result.OneCVSet;
import mdr.result.SavedModels;
import mdr.result.Suite;

import util.NewIt;

/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */

public class LinearMergeSearch extends AbstractMergeSearch {

	protected int[] SNPIndex;
	protected Combination model;
	protected int count;

	public LinearMergeSearch(DataFile dr, Subdivision sd) {
		super(dr, sd);
	}

	public void search(int or, ArrayList<String> modelspace) {
		order = or;
		bestKFold = new BestKFoldCVResult(or, subdivision.getInterval());
		count = 0;
		Arrays.fill(currBestStats, 0);

		savedModels = new SavedModels();
		ArrayList<DataFile.Subject> sample = data.getSample();
		for (String modelName : modelspace) {
			for (Combination testingModel : cvTestingSet) {
				testingModel.clear();
			}
			model = new Combination();
			SNPIndex = ToolKit.StringToIntArray(modelName);
			// mergeSearch(sample, null, 0);
			linearSearch(sample);
			calculate(modelName);
			count++;
		}
		summary();
	}

	public void linearSearch(ArrayList<DataFile.Subject> subjects) {
		for (DataFile.Subject sub : subjects) {
			String geno = sub.getGenotype(SNPIndex);
			if (geno.contains(Parameter.missing_genotype)) {
				continue;
			} else {
				Suite subset = model.get(geno);
				if (subset == null) {
					subset = new Suite();
					model.put(geno, subset);
				}
				subset.add(sub);
			}
		}

		for (Entry<String, Suite> entry : model.entrySet()) {
			String geno = entry.getKey();
			Suite s = entry.getValue();
			s.summarize();
			assignKFold(geno, s.getSubjects());
		}
	}

	public void mergeSearch(ArrayList<DataFile.Subject> subjects, String combination, int idxMarker) {
		if (idxMarker < SNPIndex.length) {
			HashMap<String, ArrayList<DataFile.Subject>> subsets = NewIt.newHashMap();
			for (DataFile.Subject sub : subjects) {
				String m = Byte.toString(sub.getGenotype(SNPIndex[idxMarker]));
				ArrayList<DataFile.Subject> subset = subsets.get(m);
				if (subset == null) {
					subset = new ArrayList<DataFile.Subject>();
					subsets.put(m, subset);
				}
				subset.add(sub);
			}
			for (String key : subsets.keySet()) {
				StringBuilder com = new StringBuilder();
				if (combination != null) {
					com.append(combination);
					com.append(MDRConstant.seperator);
				}
				com.append(key);
				mergeSearch(subsets.get(key), com.toString(), idxMarker + 1);
			}
		} // if we've processed all attributes
		else {
			run(combination, subjects);
		}
	}

	protected void run(String com, ArrayList<DataFile.Subject> subsample) {
		Suite s = new Suite(subsample);
		s.summarize();
		model.put(com, s);
		assignKFold(com, subsample);
	}

	protected void assignKFold(String key, ArrayList<DataFile.Subject> subsample) {
		// assign each individual to one's testing set
		for (DataFile.Subject sub : subsample) {
			Integer ID = sub.getIntegerID();
			int d = dataPartitionMap.get(ID).intValue();
			Combination suiteMap = cvTestingSet.get(d);
			Suite S = suiteMap.get(key);
			if (S == null) {
				S = new Suite();
				suiteMap.put(key, S);
			}
			S.add(sub);
		}

		// summarize the testing set
		for (Combination testingModels : cvTestingSet) {
			if (testingModels.containsKey(key)) {
				Suite testingSuite = testingModels.get(key);
				testingSuite.summarize();
			}
		}
	}

	public void calculate(String currModel) {
		for (int i = 0; i < subdivision.getInterval(); i++) {
			OneCVSet cvSet = new OneCVSet(i, currModel);
			Combination testingModels = cvTestingSet.get(i);
			int tr_status = 0;
			Cell trCell;
			Cell tCell;
			int idx = 0;
			for (String cellKey : model.keySet()) {
				Suite fullSuite = model.get(cellKey);
				int fullposSubs = fullSuite.getPositiveSubjects();
				int fullnegSubs = fullSuite.getNegativeSubjects();
				double fullposScr = fullSuite.getPositiveScore();
				double fullnegScr = fullSuite.getNegativeScore();
				if (testingModels.containsKey(cellKey)) {
					Suite testingSuite = testingModels.get(cellKey);
					int posSubs = testingSuite.getPositiveSubjects();
					int negSubs = testingSuite.getNegativeSubjects();
					double posScr = testingSuite.getPositiveScore();
					double negScr = testingSuite.getNegativeScore();
					int trposSubs = fullposSubs - posSubs;
					int trnegSubs = fullnegSubs - negSubs;
					double trposScr = fullposScr - posScr;
					double trnegScr = fullnegScr - negScr;
					tr_status = ToolKit.Acertainment(trposScr, trnegScr);
					trCell = new Cell(trposSubs, trnegSubs, trposScr, trnegScr, tr_status);
					tCell = new Cell(posSubs, negSubs, posScr, negScr, tr_status);
				} else {
					tr_status = ToolKit.Acertainment(fullposScr, fullnegScr);
					trCell = new Cell(fullposSubs, fullnegSubs, fullposScr, fullnegScr, tr_status);
					tCell = new Cell(0, 0, 0, 0, -1);
				}
				cvSet.addTrainingModel(cellKey, trCell);
				cvSet.addTestingModel(cellKey, tCell);
				idx++;
			}
			double trAccu = 0;
			trAccu = ToolKit.BalancedAccuracy(cvSet.getTrainingSubdivision());
			cvSet.setStatistic(MDRConstant.TrainingBalancedAccuIdx, trAccu);
			System.out.println(currModel + " " + trAccu);
			if (count == 0) {
				currBestStats[i] = trAccu;
				bestKFold.add(cvSet);
				savedModels.save(currModel, model, null);
			} else if (currBestStats[i] < trAccu) {
				currBestStats[i] = trAccu;
				String deleteModel = bestKFold.getModelAtCV(i);
				bestKFold.set(i, cvSet);
				savedModels.save(currModel, model, deleteModel);
			}
		}
	}

	private void summary() {
		for (int j = 0; j < subdivision.getInterval(); j++) {
			OneCVSet cvSet = bestKFold.get(j);
			double testingAccu = 0;
			testingAccu = ToolKit.BalancedAccuracy(cvSet.getTestingSubdivision());
			cvSet.setStatistic(MDRConstant.TestingBalancedAccuIdx, testingAccu);
		}
		bestKFold.summarise();
		stats = bestKFold.getKFoldStats();
	}

	public double[] calculateSingleBest(String modelName) {

		double[] mean = new double[MDRConstant.NumStats];
		for (int j = 0; j < subdivision.getInterval(); j++) {
			OneCVSet cvSet = new OneCVSet(j, modelName);
			Combination testingModels = (Combination) cvTestingSet.get(j);
			int tr_status;
			int t_status;
			Cell trCell;
			Cell tCell;
			Set<String> cellKeys = model.keySet();
			double[] trStatus = new double[cellKeys.size()];
			double[] tStatus = new double[cellKeys.size()];
			int idx = 0;
			for (String cellKey : cellKeys) {
				Suite fullSuite = (Suite) model.get(cellKey);
				int fullposSubs = fullSuite.getPositiveSubjects();
				int fullnegSubs = fullSuite.getNegativeSubjects();
				double fullposScr = fullSuite.getPositiveScore();
				double fullnegScr = fullSuite.getNegativeScore();
				if (testingModels.containsKey(cellKey)) {
					Suite testingSuite = (Suite) testingModels.get(cellKey);
					int posSubs = testingSuite.getPositiveSubjects();
					int negSubs = testingSuite.getNegativeSubjects();
					double posScr = testingSuite.getPositiveScore();
					double negScr = testingSuite.getNegativeScore();
					int trposSubs = fullposSubs - posSubs;
					int trnegSubs = fullnegSubs - negSubs;
					double trposScr = fullposScr - posScr;
					double trnegScr = fullnegScr - negScr;
					tr_status = ToolKit.Acertainment(trposScr, trnegScr);
					t_status = ToolKit.Acertainment(posScr, negScr);
					trCell = new Cell(trposSubs, trnegSubs, trposScr, trnegScr, tr_status);
					tCell = new Cell(posSubs, negSubs, posScr, negScr, t_status);
				} else {
					tr_status = ToolKit.Acertainment(fullposScr, fullnegScr);
					t_status = 1 - tr_status;
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

		return mean;
	}

	public double[] singleBest(String com) {
		for (int i = 0; i < cvTestingSet.size(); i++) {
			Combination testingModel = (Combination) cvTestingSet.get(i);
			testingModel.clear();
		}
		model = new Combination();
		SNPIndex = ToolKit.StringToIntArray(com);
		mergeSearch(data.getSample(), null, 0);
		return calculateSingleBest(com);
	}

	@Override
	public void summarise() {
		// TODO Auto-generated method stub

	}

	public String toString() {
		return bestKFold.getBestModel();
	}
}
