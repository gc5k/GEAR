package family.mdr.selector;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Random;

import family.mdr.arsenal.MDRStatistic;

import util.NewIt;

public class TopN implements Iterator<String> {
	
	private ArrayList<String> OrderKeyList = null;
	private HashMap<String, MDRStatistic> ModelStats = null;
	private String MaxKey = null;
	private String MinKey = null;
	private boolean flag = false;
	private int count = 0;
	private int N = 5;
	private int index = 0;

	public TopN(int n) {

		N = n;
		OrderKeyList = NewIt.newArrayList();
		OrderKeyList.ensureCapacity(N);
		ModelStats = NewIt.newHashMap();

	}

	public void add(String model, MDRStatistic stat) {

		if (!flag) {

			count++;
			ModelStats.put(model, stat);
			OrderKeyList.add(model);

			order(OrderKeyList.size() - 1);

			if (count == N)
				flag = true;

		} else {

			if (stat.compareTo(ModelStats.get(OrderKeyList.get(N - 1))) > 0) {
				
				ModelStats.remove(OrderKeyList.get(N - 1));
				ModelStats.put(model, stat);

				OrderKeyList.set(N - 1, model);
				order(N - 1);

			}

		}
		MinKey = OrderKeyList.get(count - 1);
		MaxKey = OrderKeyList.get(0);

	}

	private void order(int n) {
		for (int i = n; i >= 1; i--) {
			if (ModelStats.get(OrderKeyList.get(i)).compareTo(ModelStats.get(OrderKeyList.get(i - 1))) > 0) {
				String tmp = OrderKeyList.get(i);
				OrderKeyList.set(i, OrderKeyList.get(i - 1));
				OrderKeyList.set(i - 1, tmp);
			} else {
				break;
			}
		}
	}

	public MDRStatistic getMaxStat() {
		return ModelStats.get(MaxKey);
	}

	public MDRStatistic getMinStat() {
		return ModelStats.get(MinKey);
	}

	public String getMaxKey() {
		return MaxKey;
	}
	
	public String getMinKey() {
		return MinKey;
	}

	public MDRStatistic getMDRStatistic(String key) {
		if(ModelStats.containsKey(key)) {
			return ModelStats.get(key);
		} else {
			return null;
		}
	}

	public HashMap<String, MDRStatistic> getResult() {
		return ModelStats;
	}

	public int getKeyListLength() {
		return OrderKeyList.size();
	}
	
	@Override
	public boolean hasNext() {
		if(index == OrderKeyList.size()) {
			return false;
		} else {
			return true;
		}
	}

	@Override
	public String next() {
		if (OrderKeyList == null) {
			return null;
		} else {
			return OrderKeyList.get(index++);
		}
	}

	@Override
	public void remove() {
		throw new UnsupportedOperationException();
	}

	public static void main(String[] args) {
		String[] k = { "0", "1", "2", "3", "4", "5", "6", "7", "8", "9" };
		double[] d = new double[k.length];
		Random rand = new Random();
		TopN topN = new TopN(5);
		for (int i = 0; i < k.length; i++) {
			d[i] = rand.nextDouble();
			MDRStatistic m = new MDRStatistic(d[i], rand.nextDouble());
			topN.add(k[i], m);
		}
		System.out.println();
	}

}
