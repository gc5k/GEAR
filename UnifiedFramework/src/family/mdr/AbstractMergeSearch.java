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
	protected double Threshold;

	protected boolean mute = true;
	protected Random rnd = new Random(Parameter.seed);

	public AbstractMergeSearch(int c, ArrayList<PersonIndex> dr, MapFile mf, ModelGenerator mg, int n, boolean m) {
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
	
	abstract protected void linearSearch();

	protected void cleanupTestingSet() {
		for (Combination testingModel : cvTestingSet) {
			testingModel.clear();
		}
	}

	public void setMute(boolean flag) {
		mute = flag;
	}

	public HashMap<String, MDRStatistic> getMDRResult() {
		HashMap<String, MDRStatistic> m = NewIt.newHashMap();
		m.put(bestModel, bestStat);
		return m;
	}

	public double[] getModelStats() {
		return bestStat.getStats();
	}
	
	public abstract void calculateSingleBest(String modelName);

	public abstract void search(int or, int N);

	public abstract String toString();
	
	public abstract void kernal(String modelName);

}
