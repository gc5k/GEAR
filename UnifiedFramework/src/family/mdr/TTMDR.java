package family.mdr;

import java.util.ArrayList;
import java.util.Set;
import java.util.Map.Entry;

import statistics.FisherExactTest.mdrExactTest.MDRTestingExactTest;
import statistics.FisherExactTest.mdrExactTest.MDRTrainingExactTest;

import admixture.parameter.Parameter;

import family.mdr.arsenal.MDRConstant;
import family.mdr.arsenal.ModelGenerator;
import family.mdr.arsenal.ToolKit;
import family.mdr.data.PersonIndex;
import family.mdr.result.Cell;
import family.mdr.result.Combination;
import family.mdr.result.MDRStatistic;
import family.mdr.result.OneCVSet;
import family.mdr.result.Suite;
import family.pedigree.file.MapFile;

public class TTMDR extends AbstractMergeSearch {

	public TTMDR(int c, ArrayList<PersonIndex> dr, MapFile mf, ModelGenerator mg, int n, boolean m) {
		super(c, dr, mf, mg, n, m);

	}

	@Override
	public void search(int or, int N) {

		bestStat = new MDRStatistic();
		order = or;
		cg.revup(or);
		count = 0;
		topN = N;

		int count = 0;

		for (; cg.hasNext();) {
			String m = cg.next();
			if (rnd.nextDouble() > Parameter.thin)
				continue;
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

	@Override
	public void kernal(String modelName) {

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

	@Override
	protected void linearSearch() {
		model = new Combination();
		mdrStat = new MDRStatistic();
		double Tp = 0;
		double Tn = 0;
		double T = 1;
		int TrN = 0;

		for (PersonIndex sub : data) {
			String geno = sub.getGenotype(SNPIndex);
			if (geno.contains(MDRConstant.missingGenotype)) {
				continue;
			} else {

				double s = sub.getScore();

				if (sub.getGroup() == 0) {
					if (s > 0) {
						Tp += s;
					} else {
						Tn += s;
					}
					TrN++;
				}

				Suite subset = model.get(geno);
				if (subset == null) {
					subset = new Suite();
					model.put(geno, subset);
				}
				subset.add(sub);

				int d = sub.getGroup();
				Combination suiteMap = cvTestingSet.get(d);
				Suite S = suiteMap.get(geno);
				if (S == null) {
					S = new Suite();
					suiteMap.put(geno, S);
				}
				S.add(sub);
			}
		}
		try {
			T = -1 * Tp / Tn;
		} catch (Exception E) {
			System.err.println("Denominator is zero.");
		}
		Suite.setThreshold(T);

		for (Entry<String, Suite> entry : model.entrySet()) {
			String geno = entry.getKey();
			int c = 0;
			for (Combination testingModels : cvTestingSet) {
				if (testingModels.containsKey(geno)) {
					Suite testingSuite = testingModels.get(geno);
					testingSuite.summarize();
					if (c == 0) {
						int g = Suite.Ascertainment(testingSuite.getPositiveScore(), testingSuite.getNegativeScore());
						testingSuite.setStatus(g);
					}
					c++;
				}
			}
		}
	}

	@Override
	public void calculateSingleBest(String modelName) {

		OneCVSet cvSet = new OneCVSet(0, modelName);
		Combination trainingModels = cvTestingSet.get(0);
		Combination testingModels = cvTestingSet.get(1);

		Set<String> cellKeys = trainingModels.keySet();

		for (String cellKey : cellKeys) {
			int tr_status = -1;
			Cell trCell = null;
			Cell tCell = null;

			Suite trainingSuite = trainingModels.get(cellKey);
			int pos_Subs = trainingSuite.getPositiveSubjects();
			int neg_Subs = trainingSuite.getNegativeSubjects();
			double pos_Scr = trainingSuite.getPositiveScore();
			double neg_Scr = trainingSuite.getNegativeScore();
			tr_status = trainingSuite.getStatus();
			trCell = new Cell(pos_Subs, neg_Subs, pos_Scr, neg_Scr, tr_status);
			cvSet.addTrainingModel(cellKey, trCell);
			if (testingModels.containsKey(cellKey)) {
				Suite testingSuite = testingModels.get(cellKey);
				int tpos_Subs = testingSuite.getPositiveSubjects();
				int tneg_Subs = testingSuite.getNegativeSubjects();
				double tpos_Scr = testingSuite.getPositiveScore();
				double tneg_Scr = testingSuite.getNegativeScore();
				tCell = new Cell(tpos_Subs, tneg_Subs, tpos_Scr, tneg_Scr, tr_status);
				cvSet.addTestingModel(cellKey, tCell);
			}
		}

		cellKeys = testingModels.keySet();
		for (String cellKey : cellKeys) {
			if (!trainingModels.containsKey(cellKey)) {
				Suite testingSuite = testingModels.get(cellKey);
				int pos_Subs = testingSuite.getPositiveSubjects();
				int neg_Subs = testingSuite.getNegativeSubjects();
				double pos_Scr = testingSuite.getPositiveScore();
				double neg_Scr = testingSuite.getNegativeScore();
				Cell tCell = new Cell(pos_Subs, neg_Subs, pos_Scr, neg_Scr, -1);
				cvSet.addTestingModel(cellKey, tCell);
			}
		}

		MDRTrainingExactTest mdrTrET = new MDRTrainingExactTest(cvSet.getTrainingSubdivision());
		MDRTestingExactTest mdrTET = new MDRTestingExactTest(cvSet.getTestingSubdivision());

		double trP = mdrTrET.getOneTailP();
		double tP = mdrTET.getOneTailP();
		System.err.println(modelName + " " + trP + " " + mdrTrET.getOneTailP() + ", " + tP + " " + mdrTET.getOneTailP());

		mdrStat.setTrainingPValue(trP);
		mdrStat.setTestingBalancedAccuracy(tP);

	}

	@Override
	public String toString() {
		// TODO Auto-generated method stub
		return null;
	}

}
