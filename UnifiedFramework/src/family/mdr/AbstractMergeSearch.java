package family.mdr;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;

import family.pedigree.file.MapFile;

import admixture.parameter.Parameter;

import mdr.algorithm.Subdivision;
import mdr.data.DataFile;
import mdr.result.Combination;
import mdr.result.SavedModels;
import mdr.result.Suite;

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

	protected HashMap<Integer, Integer> dataPartitionMap;
	protected DataFile data;
	protected SavedModels savedModels;
	protected Subdivision subdivision;
	protected MapFile mapData;

	protected int[] SNPIndex;
	protected Combination model;
	protected int count;

	protected int topN = 5;

	public AbstractMergeSearch(DataFile dr, Subdivision sd, MapFile mf, int[] in_snp, int[] ex_snp, int n) {
		data = dr;
		subdivision = sd;
		mapData = mf;
		cg = new CombinationGenerator(mf.getMarkerNumber(), in_snp, ex_snp);
		for (int i = 0; i < subdivision.getInterval(); i++) {
			Combination testingMap = new Combination();
			cvTestingSet.add(testingMap);
		}
		dataPartitionMap = subdivision.getDivision();
		topN = n;
	}

	protected void linearSearch(ArrayList<DataFile.Subject> subjects) {
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

	protected void cleanupTestingSet() {
		for (Combination testingModel : cvTestingSet) {
			testingModel.clear();
		}
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

	public abstract HashMap<String, MDRStatistic> getMDRResult();

	public abstract double[] getModelStats();

	public abstract void search(int or, int N);

	public abstract String toString();

}
