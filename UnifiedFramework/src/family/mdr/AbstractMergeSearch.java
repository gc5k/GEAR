package family.mdr;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;
import java.util.Map.Entry;

import admixture.parameter.Parameter;

import family.mdr.arsenal.MDRConstant;
import family.mdr.arsenal.ModelGenerator;
import family.mdr.data.PersonIndex;
import family.mdr.result.Combination;
import family.mdr.result.MDRStatistic;
import family.mdr.result.Suite;
import family.pedigree.file.MapFile;

import statistics.FisherExactTest.mdrExactTest.MDRTruncatedExactTest;
import util.NewIt;

/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */

public abstract class AbstractMergeSearch {
	protected int order;
	protected ModelGenerator cg;
	protected ArrayList<Combination> cvTestingSet = NewIt.newArrayList();
	// always keeps the K testing models of the current combination

	protected ArrayList<PersonIndex> data;
	protected MapFile mapData;

	protected int cv;
	protected int[] SNPIndex;
	protected Combination model;
	protected int count;

	protected int topN = 5;
	protected MDRStatistic mdrStat;
	protected String bestModel;
	protected MDRStatistic bestStat;

	protected boolean mute = true;
	protected Random rnd = new Random(Parameter.seed);
	
	public AbstractMergeSearch(int c, ArrayList<PersonIndex> dr, MapFile mf,
			ModelGenerator mg, int n, boolean m) {
		cv = c;
		data = dr;
		mapData = mf;
		cg = mg;

		for (int i = 0; i < cv; i++) {
			Combination testingMap = new Combination();
			cvTestingSet.add(testingMap);
		}
		topN = n;
		mute = false;
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

		double mean = (Tp + Tn)/N;
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
			if(group == 1) {
				nP += s.getPositiveSubjects() + s.getNegativeSubjects();
				mP += s.getMeanScore() * ( s.getPositiveSubjects() + s.getNegativeSubjects());
			} else if (group == 0){
				nN += s.getNegativeSubjects() + s.getPositiveSubjects();
				mN += s.getMeanScore() * ( s.getPositiveSubjects() + s.getNegativeSubjects());
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
		MDRTruncatedExactTest et = new MDRTruncatedExactTest(model);
		mdrStat.setTruncatedFisherOneTailP(et.getOneTailP());
		mdrStat.setTruncatedFisherTwoTailP(et.getTwoTailP());
//		System.err.println("Exact: " + et.getOneTailP());
	}

	protected void cleanupTestingSet() {
		for (Combination testingModel : cvTestingSet) {
			testingModel.clear();
		}
	}

	public void setMute(boolean flag) {
		mute = flag;
	}

	public abstract HashMap<String, MDRStatistic> getMDRResult();

	public abstract double[] getModelStats();

	public abstract void search(int or, int N);

	public abstract String toString();

}
