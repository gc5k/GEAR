package family.mdr;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Set;
import java.util.Map.Entry;

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
import family.pedigree.design.hierarchy.ChenInterface;
import family.pedigree.file.MapFile;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;

import test.Test;

public class HeteroCombinationSearchP extends AbstractMergeSearch {

	private double seT;
	private double meanT;
	private double TStat;
	private double pT = Double.NaN;
//	private double[] roundBest = null;
	private NormalDistribution nd = new NormalDistributionImpl();
	private ChenInterface chen;

	public static class Builder {

		private int cv;
		private ArrayList<PersonIndex> dr;
		private MapFile mf;
		private ModelGenerator mg;
		private int N = 1;
		private boolean mute = false;
		private ChenInterface chen;

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

		public Builder chen(ChenInterface chen) {
			this.chen = chen;
			return this;
		}

		public HeteroCombinationSearchP build() {
			return new HeteroCombinationSearchP(this);
		}
	}

	public HeteroCombinationSearchP(Builder builder) {

		super(builder.cv, builder.dr, builder.mf, builder.mg, builder.N, builder.mute);
		this.chen = builder.chen;
//		if (Parameter.permFlag) {
//			roundBest = new double[Parameter.perm];
//			Arrays.fill(roundBest, 0);
//		}

	}

	@Override
	public void calculateSingleBest(String modelName) {

		double[] mean = new double[MDRConstant.NumStats];

		for (Entry<String, Suite> entry : model.entrySet()) {
			String geno = entry.getKey();
			Suite s = entry.getValue();
			s.summarize();
			int group = Suite.Ascertainment(s.getPositiveScore(), s.getNegativeScore());
			if (group == 1) {
				s.setStatus(group);
			} else if (group == 0) {
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

		mdrStat.setTrainingBalancedAccuracy(mean[MDRConstant.TrainingBalancedAccuIdx]);
		mdrStat.setTestingBalancedAccuracy(mean[MDRConstant.TestingBalancedAccuIdx]);
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

	private void pkernal(String modelName) {
		pcalculateSingleBest(modelName);
	}

	@Override
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
	}

	@Override
	public void search(int or, int N) {
		order = or;
		PrintHeader();
		bestStat = new MDRStatistic();

		cg.revup(or);
		count = 0;
		topN = N;

		long t0 = System.currentTimeMillis();
		int c = 0;
		System.err.println("number of tests: " + cg.getTotal());
		Test.LOG.append("number of tests: " + cg.getTotal());
		Test.LOG.append("\n");
		for (; cg.hasNext();) {
			String m = cg.next();
			c++;
			if (c == Parameter.testUnit) {
				testdrive(t0);
				if(Parameter.testdrive) {
					System.exit(1);
				}
			}
			if (rnd.nextDouble() > Parameter.thin)
				continue;
			chen.RecoverScore();
			chen.recoverSeed();
			kernal(m);
			if (Parameter.permFlag)
				pkernal(m);
			count++;
			if (!mute) {
				boolean flag = true;
				if (Parameter.epFlag) {
					if (Parameter.permFlag && mdrStat.getTestingBalancedPT() > Parameter.ep) {
						flag = false;
					}
				}
				if (flag && Parameter.trainingFlag) {
					if (mdrStat.getTrainingBalancedAccuracy() < Parameter.threshold_training) {
						flag = false;
					}
				}
				if (flag && Parameter.testingFlag) {
					if (mdrStat.getTestingBalancedAccuracy() < Parameter.threshold_testing) {
						flag = false;
					}
				}
				if (flag && Parameter.vcFlag) {
					if (mdrStat.getVc() < Parameter.vc) {
						flag = false;
					}
				}

				if (!flag) {
					continue;
				}
				int[] idx = ToolKit.StringToIntArray(m);
				int len = idx.length;
				for (int j = 0; j < len; j++) {
					System.out.print(mapData.getSNP(idx[j]).getName());
					if (j != len - 1) {
						System.out.print("-");
					} else {
						System.out.print(", ");
					}
				}
				System.out.print(mdrStat);
				if (Parameter.verboseFlag) {
					System.out.print(", ");
					for (int j = 0; j < len; j++) {
						System.out.print(mapData.getSNP(idx[j]));
						if (j != idx.length - 1)
							System.out.print(", ");
					}
					System.out.print(", ");
					System.out.print(model.printModel(idx, mapData));
				}
				System.out.println();
			}
		}
//		calculateThreshold();
	}

	private void pcalculateSingleBest(String modelName) {

		if (!Parameter.permFlag) {
			return;
		}

		double tSq = 0;
		double tSum = 0;

		double[] permuvalue = new double[Parameter.perm];
		for (int i = 0; i < Parameter.perm; i++) {

			chen.getPermutedScore(Parameter.permu_scheme);
			for (Entry<String, Suite> entry : model.entrySet()) {
				String geno = entry.getKey();
				Suite s = entry.getValue();
				s.summarize();
				int group = Suite.Ascertainment(s.getPositiveScore(), s.getNegativeScore());
				if (group == 1) {
					s.setStatus(group);
				} else if (group == 0) {
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

			double[] Pmean = new double[MDRConstant.NumStats];

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

				double trAccu = 0;
				double tAccu = 0;
				trAccu = ToolKit.BalancedAccuracy(cvSet.getTrainingSubdivision());
				Pmean[MDRConstant.TrainingBalancedAccuIdx] += trAccu;

				tAccu = ToolKit.BalancedAccuracy(cvSet.getTestingSubdivision());
				Pmean[MDRConstant.TestingBalancedAccuIdx] += tAccu;
			}

			Pmean[MDRConstant.TrainingBalancedAccuIdx] /= cv;
			Pmean[MDRConstant.TestingBalancedAccuIdx] /= cv;
			permuvalue[i] = Pmean[MDRConstant.TestingBalancedAccuIdx];
			tSq += Pmean[MDRConstant.TestingBalancedAccuIdx] * Pmean[MDRConstant.TestingBalancedAccuIdx];
			tSum += Pmean[MDRConstant.TestingBalancedAccuIdx];

//			if (roundBest[i] < permuvalue[i]) {
//				roundBest[i] = permuvalue[i];
//			}
		}

		meanT = tSum / Parameter.perm;
		seT = Math.sqrt((tSq - (tSum / Parameter.perm) * tSum) / Parameter.perm);
		TStat = (mdrStat.getTestingBalancedAccuracy() - meanT) / seT;
		mdrStat.setTestingBalancedAccuracyMeanT(meanT);
		mdrStat.setTestingBalancedAccuracyseT(seT);
		try {
			pT = 1 - nd.cumulativeProbability(TStat);
			mdrStat.setTestingBalancedAccuracyPT(pT);
		} catch (MathException e) {
			e.printStackTrace();
		}

	}

//	private void calculateThreshold() {
//		if (roundBest == null) {
//			return;
//		}
//		Arrays.sort(roundBest);
//		StringBuilder sb = new StringBuilder(Parameter.out);
//		sb.append(order);
//		if (Parameter.sliceFlag) {
//			sb.append(".slice" + Parameter.slice + "." + Parameter.sliceN);
//		}
//		sb.append(".thres");
//		PrintStream PW = null;
//		try {
//			PW = new PrintStream(sb.toString());
//		} catch (FileNotFoundException e) {
//			e.printStackTrace();
//		}
//		for (int i = 0; i < roundBest.length; i++) {
//			PW.println(roundBest[i]);
//		}
//		PW.close();
////		int idx_005 = (int) (roundBest.length * 0.95);
////		int idx_001 = (int) (roundBest.length * 0.99);
////		System.out.println("threshold at significance level 0.05 = " + roundBest[idx_005]);
////		System.out.println("threshold at significance level 0.01 = " + roundBest[idx_001]);
//	}

	public void PrintHeader() {
		System.out.print("model" + ", ");
		System.out.print("effective individuals, ");
		System.out.print("vc, vx, vt, F, PF(d1,d2), ");
		for (int i = 0; i < MDRConstant.NumStats; i++) {
			if (i != MDRConstant.NumStats - 1) {
				System.out.print(MDRConstant.TestStatistic[i] + "(TA), ");
				if (Parameter.permFlag) {
					System.out.print("mean TA (null dis), " + "S.E.TA (null dis), ");
				}
			} else {
				System.out.print(MDRConstant.TestStatistic[i]);
			}
		}
		if (Parameter.verboseFlag) {
			for (int i = 0; i < order; i++) {
				System.out.print(", LOCUS" + (i + 1) + ", ");
				System.out.print("CHR" + (i + 1) + ", ");
				System.out.print("POS" + (i + 1) + ", ");
				if (i < order - 1) {
					System.out.print("MAF" + (i + 1) + ", ");
				} else {
					System.out.print("MAF" + (i + 1) + ", ");
				}
			}
			System.out
					.print("classification (genotype, risk group, positive scores, positive subjects, negative score, negative subjects)");
		}
		System.out.println();
	}

	public void testdrive(long t0) {
		long t1 = System.currentTimeMillis();
		long t = t1 - t0;
		BigInteger bi = cg.getTotal();
		double b0 = bi.longValue();
		b0 /= Parameter.testUnit;
		b0 *= t;
		b0 /= 3600000;
		System.err.println("Estimated time (hour) to complete the whole analysis: " + String.format("%.4f", b0));
		Test.LOG.append("Estimated time (hour) to complete the whole analysis: " + String.format("%.4f", b0));
		Test.LOG.append("\n");
	}

	@Override
	public String toString() {
		return null;
	}

}
