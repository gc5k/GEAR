package family.mdr;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;

import family.mdr.data.MDRConstant;
import family.mdr.data.PersonIndex;
import family.mdr.result.Combination;
import family.mdr.result.SavedModels;
import family.mdr.result.Suite;
import family.pedigree.file.MapFile;

import util.NewIt;

/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */

public abstract class AbstractMergeSearch {
	protected int order;
	protected CombinationGenerator cg;
	protected ArrayList<Combination> cvTestingSet = NewIt.newArrayList();
	// always keeps the K testing models of the current combination

	protected ArrayList<PersonIndex> data;
	protected SavedModels savedModels;
	protected MapFile mapData;

	protected int cv;
	protected int[] SNPIndex;
	protected Combination model;
	protected int count;

	protected int topN = 5;
	MDRStatistic mdrstat;

	protected boolean mute = true;

	public AbstractMergeSearch(int c, ArrayList<PersonIndex> dr, MapFile mf, int[] in_snp, int[] ex_snp, int n, boolean m) {
		cv = c;
		data = dr;
		mapData = mf;
		cg = new CombinationGenerator(mf.getMarkerNumber(), in_snp, ex_snp);
		for (int i = 0; i < cv; i++) {
			Combination testingMap = new Combination();
			cvTestingSet.add(testingMap);
		}
		topN = n;
		mute = false;
	}

	protected void linearSearch() {
		for (PersonIndex sub : data) {
			String geno = sub.getGenotype(SNPIndex);
			if (geno.contains(MDRConstant.missingGenotype)) {
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

	protected void cleanupTestingSet() {
		for (Combination testingModel : cvTestingSet) {
			testingModel.clear();
		}
	}

	protected void assignKFold(String key, ArrayList<PersonIndex> subsample) {
		// assign each individual to one's testing set
		for (PersonIndex sub : subsample) {
			int d = sub.getGroup();
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

	public void setMute(boolean flag) {
		mute = flag;
	}

	public abstract HashMap<String, MDRStatistic> getMDRResult();

	public abstract double[] getModelStats();

	public abstract void search(int or, int N);

	public abstract String toString();

}
