package family.mdr;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;

import mdr.MDRConstant;
import mdr.algorithm.Subdivision;
import mdr.arsenal.ToolKit;
import mdr.data.DataFile;
import mdr.result.Cell;
import mdr.result.Combination;
import mdr.result.OneCVSet;
import mdr.result.Suite;
import family.pedigree.file.MapFile;

public class HeteroCombinationSearch extends AbstractMergeSearch {

	private double[] statistic;
	private TopN Top_N = null;
	
	public static class Builder {
		private DataFile dr;
		private Subdivision sd;
		private MapFile mf;
		int[] in_snp = null;
		int[] ex_snp = null;
		
		int N = 1;
		
		public Builder(DataFile dr, Subdivision sd, MapFile mf, int[] in_snp, int[] ex_snp) {
			this.dr = dr;
			this.sd = sd;
			this.mf = mf;
			this.in_snp = in_snp;
			this.ex_snp = ex_snp;
		}
		
		public Builder topN(int n) {
			this.N = n;
			return this;
		}

		public HeteroCombinationSearch build() {
			return new HeteroCombinationSearch(this);
		}
	}
	
	public HeteroCombinationSearch(Builder builder) {
		super(builder.dr, builder.sd, builder.mf, builder.in_snp, builder.ex_snp, builder.N);
	}

	public void search(int or, int N) {
		statistic = new double[MDRConstant.NumStats];

		order = or;
		cg.revup(or);
		count = 0;
		topN = N;

		Top_N = new TopN(topN);

		int count = 0;

		for( ; cg.hasNext(); ) {
				kernal(cg.next());
				count++;
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

		mdrstat.setTrainingBalancedAccuracy(mean[MDRConstant.TrainingBalancedAccuIdx]);
		mdrstat.setTestingBalancedAccuracy(mean[MDRConstant.TestingBalancedAccuIdx]);
		if(count<topN) {
			Top_N.add(modelName, mdrstat);
		} else {
			if(mdrstat.compareTo(Top_N.getMinStat())>0) {
				Top_N.add(modelName, mdrstat);
			}
		}
		return mean;
	}

	public double[] getStats() {
		return statistic;
	}

	public double[] getModelStats() {
		MDRStatistic m = Top_N.getMDRStatistic(Top_N.getMaxKey());
		return m.getStats();
		
	}
	public HashMap<String, MDRStatistic> getMDRResult() {
		return Top_N.getResult();
	}

	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("The top " + Top_N.getKeyListLength() + " models under " + order + " order interaction");
		sb.append(System.getProperty("line.separator"));
		sb.append("model(marker chr pos): ");
		for (int i = 0; i < MDRConstant.NumStats; i++) {
			if( i != MDRConstant.NumStats - 1) {
				sb.append(MDRConstant.TestStatistic[i] + ", ");
			} else {
				sb.append(MDRConstant.TestStatistic[i]);
			}
		}
		sb.append(System.getProperty("line.separator"));
		for (; Top_N.hasNext(); ) {
			String key = Top_N.next();
			int[] idx = ToolKit.StringToIntArray(key);
			for (int j = 0; j < idx.length; j++) {
				sb.append(mapData.getSNP(idx[j]));
				if (j != idx.length - 1)
					sb.append(", ");
			}
			sb.append(": ");
			sb.append(Top_N.getMDRStatistic(key));
			sb.append(System.getProperty("line.separator"));
		}
		sb.append("------------------------");
		return sb.toString();
	}
}
