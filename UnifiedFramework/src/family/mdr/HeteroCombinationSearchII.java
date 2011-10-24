package family.mdr;

import java.util.ArrayList;
import java.util.Set;
import java.util.Map.Entry;

import admixture.parameter.Parameter;

//import statistics.FisherExactTest.mdrExactTest.MDRTestingExactTest;
//import statistics.FisherExactTest.mdrExactTest.MDRTrainingExactTest;
//import statistics.FisherExactTest.mdrExactTest.MDRTruncatedExactTest;

import family.mdr.arsenal.MDRConstant;
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
			if (rnd.nextDouble() > Parameter.thin)
				continue;
			kernal(m);
			count++;
			if (!mute) {
				boolean flag = true;
				if (Parameter.epFlag ) {
					if (mdrStat.getTestingBalancedAccuracy() < Parameter.ep) {
						flag = false;
					}
				} else {
					if (Parameter.trainingFlag) {
						if (mdrStat.getTrainingBalancedAccuracy() < Parameter.threshold_training ) {
							flag = false;
						}
					}
					if (Parameter.testingFlag) {
						if (mdrStat.getTestingBalancedAccuracy() < Parameter.threshold_testing) {
							flag = false;
						}
					}
				}
				if (!flag ) {
					continue;
				}
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

	protected void linearSearch() {

		model = new Combination();
		mdrStat = new MDRStatistic();
		double Tp = 0;
		double Tn = 0;
		double T = 1;
		int N = 0;
		double Vt = 0;

		for (PersonIndex sub : data) {
			String geno = sub.getGenotype(SNPIndex);
			if (geno.contains(MDRConstant.missingGenotype)) {
				continue;
			} else {

				double s = sub.getScore();
				N++;
				Vt += s * s;
				if (s > 0) {
					Tp += s;
				} else {
					Tn += s;
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
			T /= -1 * Tp / Tn;
		} catch (Exception E) {
			System.err.println("Denominator is zero.");
		}
		Suite.setThreshold(T);

		double mean = (Tp + Tn) / N;
		Vt -= N * mean * mean;
		mdrStat.setVt(Vt);
		mdrStat.setN(N);

		int nP = 0;
		int nN = 0;
		double mP = 0;
		double mN = 0;

		for (Entry<String, Suite> entry : model.entrySet()) {
			String geno = entry.getKey();
			Suite s = entry.getValue();
			s.summarize();
			int group = Suite.Ascertainment(s.getPositiveScore(), s.getNegativeScore());
			if (group == 1) {
				nP += s.getPositiveSubjects() + s.getNegativeSubjects();
				mP += s.getMeanScore() * (s.getPositiveSubjects() + s.getNegativeSubjects());
				s.setStatus(group);
			} else if (group == 0) {
				nN += s.getNegativeSubjects() + s.getPositiveSubjects();
				mN += s.getMeanScore() * (s.getPositiveSubjects() + s.getNegativeSubjects());
				s.setStatus(group);
			} else {

			}
			for (Combination testingModels : cvTestingSet) {
				if (testingModels.containsKey(geno)) {
					Suite testingSuite = testingModels.get(geno);
					testingSuite.summarize();
				}
			}
		}

		double meanPos = 0;
		double meanNeg = 0;
		if (mP != 0 && mN != 0) {
			meanPos = mP / nP;
			meanNeg = mN / nN;
		} else if (mP != 0 && mN == 0) {
			meanPos = mP / nP;
		} else if (mP == 0 && mN != 0) {
			meanNeg = mN / nN;
		}
		mdrStat.setNpos(nP);
		mdrStat.setNneg(nN);
		double Vx = nP * (meanPos - mean) * (meanPos - mean) + nN * (meanNeg - mean) * (meanNeg - mean);
		mdrStat.setVx(Vx);
//		MDRTruncatedExactTest et = new MDRTruncatedExactTest(model);
//		mdrStat.setTruncatedFisherOneTailP(et.getOneTailP());
//		mdrStat.setTruncatedFisherTwoTailP(et.getTwoTailP());
//		System.err.println("Exact: " + et.getOneTailP());
	}

	public void calculateSingleBest(String modelName) {

		double[] mean = new double[MDRConstant.NumStats];
		for (int j = 0; j < cv; j++) {
			OneCVSet cvSet = new OneCVSet(j, modelName);
			Combination testingModels = cvTestingSet.get(j);

			Cell trCell;
			Cell tCell;
			Set<String> cellKeys = model.keySet();

			for (String cellKey : cellKeys) {
				int tr_status = -1;
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
					tr_status = Suite.Ascertainment(fullposScr - pos_Scr, fullnegScr - neg_Scr);
					trCell = new Cell(fullposSubs - pos_Subs, fullnegSubs - neg_Subs, fullposScr - pos_Scr, fullnegScr - neg_Scr, tr_status);
					tCell = new Cell(pos_Subs, neg_Subs, pos_Scr, neg_Scr, tr_status);
				} else {
					tr_status = Suite.Ascertainment(fullposScr, fullnegScr);
					trCell = new Cell(fullposSubs, fullnegSubs, fullposScr, fullnegScr, tr_status);
					tCell = new Cell(0, 0, 0, 0, -1);
				}
				cvSet.addTrainingModel(cellKey, trCell);
				cvSet.addTestingModel(cellKey, tCell);

			}
			
//			MDRTrainingExactTest mdrTrET = new MDRTrainingExactTest(cvSet.getTrainingSubdivision());
//			MDRTestingExactTest mdrTET = new MDRTestingExactTest(cvSet.getTestingSubdivision());

			double trAccu = 0;
			double tAccu = 0;
			trAccu = ToolKit.BalancedAccuracy(cvSet.getTrainingSubdivision());
			mean[MDRConstant.TrainingBalancedAccuIdx] += trAccu;

			tAccu = ToolKit.BalancedAccuracy(cvSet.getTestingSubdivision());
			mean[MDRConstant.TestingBalancedAccuIdx] += tAccu;
//			System.err.println(modelName + " " + trAccu + " " + mdrTrET.getOneTailP()+ ", " + tAccu + " " + mdrTET.getOneTailP());
			cvSet.setStatistic(MDRConstant.TrainingBalancedAccuIdx, trAccu);
			cvSet.setStatistic(MDRConstant.TestingBalancedAccuIdx, tAccu);
		}

		mean[MDRConstant.TrainingBalancedAccuIdx] /= cv;
		mean[MDRConstant.TestingBalancedAccuIdx] /= cv;

		mdrStat.setTrainingBalancedAccuracy(mean[MDRConstant.TrainingBalancedAccuIdx]);
		mdrStat.setTestingBalancedAccuracy(mean[MDRConstant.TestingBalancedAccuIdx]);
	}

	public void PrintHeader() {
		System.out.print(System.getProperty("line.separator"));
		System.out.print("model code, model(marker chr pos minor allele major allele): ");
		for (int i = 0; i < MDRConstant.NumStats; i++) {
			if (i != MDRConstant.NumStats - 1) {
				System.out.print(MDRConstant.TestStatistic[i] + ", ");
			} else {
				System.out.print(MDRConstant.TestStatistic[i]);
			}
		}
		System.out.print(", Truncated Fisher's Exact Test" + ", ");
		System.out
				.print(": classfication (genotype, High-risk or Low-risk group, positive scores, positive subjects, negative score, negative subjects)");
		System.out.println();
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
