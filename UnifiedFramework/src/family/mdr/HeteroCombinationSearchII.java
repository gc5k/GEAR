package family.mdr;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;

import admixture.parameter.Parameter;

import util.NewIt;

import family.mdr.arsenal.MDRConstant;
import family.mdr.arsenal.MDRStatistic;
import family.mdr.arsenal.ModelGenerator;
import family.mdr.arsenal.ToolKit;
import family.mdr.data.PersonIndex;
import family.mdr.result.Cell;
import family.mdr.result.Combination;
import family.mdr.result.OneCVSet;
import family.mdr.result.Suite;
import family.pedigree.file.MapFile;

public class HeteroCombinationSearchII extends AbstractMergeSearch {

	public static class Builder {

		private int cv;
		private ArrayList<PersonIndex> dr;
		private MapFile mf;
		private ModelGenerator mg;
		private int[][] wseq2 = null;
		private int N = 1;
		private boolean mute = false;

		public Builder(int cv, ArrayList<PersonIndex> dr, MapFile mf) {
			this.cv = cv;
			this.dr = dr;
			this.mf = mf;
		}
		
		public Builder ModelGenerator(ModelGenerator mg) {
			this.mg = mg;
			return this;
		}

		public Builder mute(boolean m) {
			this.mute = m;
			return this;
		}

		public HeteroCombinationSearchII build() {
			return new HeteroCombinationSearchII(this);
		}
	}

	public HeteroCombinationSearchII(Builder builder) {

		super(builder.cv, builder.dr, builder.mf, builder.mg, builder.N, builder.mute);
	}

	public void search(int or, int N) {

		bestStat = new MDRStatistic();
		order = or;
		cg.revup(or);
		count = 0;
		topN = N;

		int count = 0;

		for (; cg.hasNext();) {
			String m = cg.next();
			if (rnd.nextDouble() > Parameter.thin) continue;
			kernal(m);
			count++;
			if (!mute) {
				int[] idx = ToolKit.StringToIntArray(m);
				System.out.print(m + ", ");
				for (int j = 0; j < idx.length; j++) {
					System.out.print(mapData.getSNP(idx[j]));
					if (j != idx.length - 1)
						System.out.print(", ");
				}
				System.out.print(";\t");
				System.out.print(mdrStat);
				System.out.print(";\t");
				System.out.print(model.printModel(idx, mapData));
				System.out.print(System.getProperty("line.separator"));
			}
		}
	}

	private void kernal(String modelName) {

		cleanupTestingSet();

		SNPIndex = ToolKit.StringToIntArray(modelName);
		linearSearch();

		calculateSingleBest(modelName);
		if (mdrStat.compareTo(bestStat) > 0) {
			bestModel = modelName;
			bestStat = mdrStat;
		}
		count++;
	}

	public void calculateSingleBest(String modelName) {

		double[] mean = new double[MDRConstant.NumStats];
		for (int j = 0; j < cv; j++) {
			OneCVSet cvSet = new OneCVSet(j, modelName);
			Combination testingModels = cvTestingSet.get(j);

			Cell trCell;
			Cell tCell;
			Set<String> cellKeys = model.keySet();
			double[] trStatus = new double[cellKeys.size()];
			double[] tStatus = new double[cellKeys.size()];
			int idx = 0;
			for (String cellKey : cellKeys) {
				int tr_status = -1;
				int t_status = -1;
				Suite fullSuite = model.get(cellKey);
				int fullposSubs = fullSuite.getPositiveSubjects();
				int fullnegSubs = fullSuite.getNegativeSubjects();
				double fullposScr = fullSuite.getPositiveScore();
				double fullnegScr = fullSuite.getNegativeScore();
				if (testingModels.containsKey(cellKey)) {
					Suite testingSuite = testingModels.get(cellKey);
					int pos_Subs = testingSuite.getPositiveSubjects();
					int neg_Subs = testingSuite.getNegativeSubjects();
					double pos_Scr = testingSuite.getPositiveScore();
					double neg_Scr = testingSuite.getNegativeScore();
					tr_status = Suite.Ascertainment(fullposScr - pos_Scr,
							fullnegScr - neg_Scr);
					t_status = Suite.Ascertainment(pos_Scr, neg_Scr);
					trCell = new Cell(fullposSubs - pos_Subs, fullnegSubs
							- neg_Subs, fullposScr - pos_Scr, fullnegScr
							- neg_Scr, tr_status);
					tCell = new Cell(pos_Subs, neg_Subs, pos_Scr, neg_Scr,
							tr_status);
				} else {
					tr_status = Suite.Ascertainment(fullposScr, fullnegScr);
					trCell = new Cell(fullposSubs, fullnegSubs, fullposScr,
							fullnegScr, tr_status);
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

		mean[MDRConstant.TrainingBalancedAccuIdx] /= cv;
		mean[MDRConstant.TestingBalancedAccuIdx] /= cv;


		mdrStat
				.setTrainingBalancedAccuracy(mean[MDRConstant.TrainingBalancedAccuIdx]);
		mdrStat
				.setTestingBalancedAccuracy(mean[MDRConstant.TestingBalancedAccuIdx]);
	}

	public double[] getModelStats() {
		return bestStat.getStats();
	}

	public HashMap<String, MDRStatistic> getMDRResult() {
		HashMap<String, MDRStatistic> m = NewIt.newHashMap();
		m.put(bestModel, bestStat);
		return m;
	}

	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("The top model under " + order + " order interaction");
		sb.append(System.getProperty("line.separator"));
		sb.append("model(marker chr pos): ");
		for (int i = 0; i < MDRConstant.NumStats; i++) {
			if (i != MDRConstant.NumStats - 1) {
				sb.append(MDRConstant.TestStatistic[i] + ", ");
			} else {
				sb.append(MDRConstant.TestStatistic[i]);
			}
		}
		sb.append(System.getProperty("line.separator"));

		String key = bestModel;
		int[] idx = ToolKit.StringToIntArray(key);
		for (int j = 0; j < idx.length; j++) {
			sb.append(mapData.getSNP(idx[j]));
			if (j != idx.length - 1)
				sb.append(", ");
		}
		sb.append(": ");
		sb.append(mdrStat);
		sb.append(System.getProperty("line.separator"));
		sb.append("------------------------");
		return sb.toString();
	}
}
