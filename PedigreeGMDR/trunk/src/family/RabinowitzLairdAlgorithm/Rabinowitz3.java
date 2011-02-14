package family.RabinowitzLairdAlgorithm;

import java.util.Iterator;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import util.NewIt;

/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class Rabinowitz3 extends AbstractGenoDistribution {
	public Rabinowitz3(TreeMap<String, Integer> children) {
		super(children);
		countChildrenAllele(childrenGenoMap);
		countAllele(childrenGenoMap);
	}

	protected void genotypeParents() {
		// TODO Auto-generated method stub
	}

	public String[] getNontransmitted(final String transmitted) {
		return null;
	}

	//
	// private void print(String[] control) {
	// for (int i = 0; i < control.length; i++) {
	// System.out.print(control[i] + "\t");
	// }
	// System.out.println();
	// }

	public String[] getNontransmitted() {
		String[] control = new String[getChildrenNum()];
		TreeMap<String, Integer> controlMap = NewIt.newTreeMap();
		if (childrenGenoMap.size() == 1) {// situation 1
			// System.err.println("Rabinowitz Table3 s1");

			String[] genopool = new String[1];
			genopool[0] = childrenGenoMap.firstKey();
			double[] freq = new double[1];
			freq[0] = 1.0;
			Produce(control, controlMap, genopool, freq);
		} else if (childrenGenoMap.size() == 2) {
			if (numHomozygous(childrenGenoMap) == 0) {
				if (numChildrenAllele() == 3) {// 5
					// System.err.println("Rabinowitz Table3 s5");

					String[] genopool = new String[2];
					genopool[0] = childrenGenoMap.firstKey();
					genopool[1] = childrenGenoMap.lastKey();
					double[] freq = new double[2];
					freq[0] = 0.5;
					freq[1] = 1;
					do {
						controlMap.clear();
						Produce(control, controlMap, genopool, freq);
					} while (!Criteria3_5(controlMap));
				} else {// 7
					// System.err.println("Rabinowitz Table3 s7");

					String[] genopool = new String[2];
					genopool[0] = childrenGenoMap.firstKey();
					genopool[1] = childrenGenoMap.lastKey();
					double[] freq = new double[2];
					freq[0] = 0.5;
					freq[1] = 1;
					do {
						controlMap.clear();
						Produce(control, controlMap, genopool, freq);
					} while (!Criteria3_7(controlMap));
				}
			}
			if (numHomozygous(childrenGenoMap) == 1) {
				if (numChildrenAllele() == 2) {// 2
					// System.err.println("Rabinowitz Table3 s2");

					shuffle(control);
				} else {// 6-1
					// System.err.println("Rabinowitz Table3 s6-1");

					String[] genopool = new String[4];
					genopool[0] = childrenGenoMap.firstKey();
					genopool[1] = childrenGenoMap.lastKey();

					char allele[][] = new char[2][2];
					if (genopool[0].charAt(0) <= genopool[1].charAt(0)) {
						allele[0][0] = genopool[0].charAt(0);
						allele[0][1] = genopool[1].charAt(0);
					} else {
						allele[0][0] = genopool[1].charAt(0);
						allele[0][1] = genopool[0].charAt(0);
					}
					if (genopool[0].charAt(1) <= genopool[1].charAt(1)) {
						allele[1][0] = genopool[0].charAt(1);
						allele[1][1] = genopool[1].charAt(1);
					} else {
						allele[1][0] = genopool[1].charAt(1);
						allele[1][1] = genopool[0].charAt(1);
					}
					genopool[2] = new String(allele[0]);
					genopool[3] = new String(allele[1]);
					double[] freq = new double[4];
					freq[0] = 0.25;
					freq[1] = 0.5;
					freq[2] = 0.75;
					freq[3] = 1;
					do {
						controlMap.clear();
						Produce(control, controlMap, genopool, freq);
					} while (!Criteria3_6(controlMap));
				}
			}
			if (numHomozygous(childrenGenoMap) == 2) {// 3-1
				// System.err.println("Rabinowitz Table3 s3-1");

				String[] genopool = new String[3];
				genopool[0] = childrenGenoMap.firstKey();
				genopool[1] = childrenGenoMap.lastKey();
				char allele[] = new char[2];
				if (genopool[0].charAt(0) <= genopool[1].charAt(0)) {
					allele[0] = genopool[0].charAt(0);
					allele[1] = genopool[1].charAt(0);
				} else {
					allele[0] = genopool[1].charAt(0);
					allele[1] = genopool[0].charAt(0);
				}
				String geno = new String(allele);
				genopool[2] = geno;
				double[] freq = new double[3];
				freq[0] = 0.25;
				freq[1] = 0.5;
				freq[2] = 1;
				do {
					controlMap.clear();
					Produce(control, controlMap, genopool, freq);
				} while (!Criteria3_3(controlMap));
			}
		} else if (childrenGenoMap.size() == 3) {
			if (numHomozygous(childrenGenoMap) == 0) {
				if (numChildrenAllele() == 3) {// 4
					// System.err.println("Rabinowitz Table3 s4");

					String[] genopool = new String[3];

					int index = 0;
					for (String g : childrenGenoMap.keySet()) {
						genopool[index] = g;
						index++;
					}
					double[] freq = new double[3];
					freq[0] = 0.33333;
					freq[1] = 0.66667;
					freq[2] = 1.0;
					do {
						controlMap.clear();
						Produce(control, controlMap, genopool, freq);
					} while (!Criteria3_4(controlMap));
				} else {// 8-1
					// System.err.println("Rabinowitz Table3 s8-1");

					String[] genopool = new String[4];
					int index = 0;
					for (String g : childrenGenoMap.keySet()) {
						genopool[index] = g;
						index++;
					}
					genopool[3] = CompatibleGenotype();
					double[] freq = new double[4];
					freq[0] = 0.25;
					freq[1] = 0.5;
					freq[2] = 0.75;
					freq[3] = 1.0;
					do {
						controlMap.clear();
						Produce(control, controlMap, genopool, freq);
					} while (!Criteria3_8(controlMap));
				}
			}
			if (numHomozygous(childrenGenoMap) == 1) {// 6-2,6-3
				// System.err.println("Rabinowitz Table3 s6-2,6-3");

				String[] genopool = new String[4];
				int index = 0;
				String[] tempgeno = new String[2];
				int ind = 0;
				for (String g : childrenGenoMap.keySet()) {
					genopool[index] = g;
					if (isHeterozygous(genopool[index])) {
						tempgeno[ind++] = genopool[index];
					}
				}
				String geno = ExtractUniqueAllele2Genotype(tempgeno[0], tempgeno[1]);
				genopool[3] = geno;

				double[] freq = new double[4];
				freq[0] = 0.25;
				freq[1] = 0.5;
				freq[2] = 0.75;
				freq[3] = 1.0;
				do {
					controlMap.clear();
					Produce(control, controlMap, genopool, freq);
				} while (!Criteria3_6(controlMap));
			}
			if (numHomozygous(childrenGenoMap) == 2) {// 3-2
				// System.err.println("Rabinowitz Table3 s3-2");

				String[] genopool = new String[3];
				int index = 0;
				for (String g:childrenGenoMap.keySet()) {
					genopool[index++] = g;
				}
				double[] freq = new double[3];
				freq[0] = 0.25;
				freq[1] = 0.5;
				freq[2] = 1.0;
				do {
					controlMap.clear();
					Produce(control, controlMap, genopool, freq);
				} while (!Criteria3_3(controlMap));
			}
		} else if (childrenGenoMap.size() == 4) {
			if (numHomozygous(childrenGenoMap) == 1) {// 6-4
				// System.err.println("Rabinowitz Table3 s6-4");

				String[] genopool = new String[4];
				int index = 0;
				for (String g:childrenGenoMap.keySet()) {
					genopool[index++] = g;
				}
				double[] freq = new double[4];
				freq[0] = 0.25;
				freq[1] = 0.5;
				freq[2] = 0.75;
				freq[3] = 1.0;
				do {
					controlMap.clear();
					Produce(control, controlMap, genopool, freq);
				} while (!Criteria3_6(controlMap));
			} else {// 8-2
				// System.err.println("Rabinowitz Table3 s8-2");

				String[] genopool = new String[4];

				int index = 0;
				for (String g:childrenGenoMap.keySet()) {
					genopool[index++] = g;
				}
				double[] freq = new double[4];
				freq[0] = 0.25;
				freq[1] = 0.5;
				freq[2] = 0.75;
				freq[3] = 1.0;
				do {
					controlMap.clear();
					Produce(control, controlMap, genopool, freq);
				} while (!Criteria3_8(controlMap));
			}
		} else {
			System.err.println("Wrecked in Rabinowitz table 3");
		}
		// print(control);
		return control;
	}

	boolean Criteria3_3(TreeMap<String, Integer> controlMap) {
		if (numHomozygous(controlMap) == 2) {
			return true;
		} else {
			return false;
		}
	}

	boolean Criteria3_4(TreeMap<String, Integer> controlMap) {
		if (controlMap.size() == 3) {
			return true;
		} else {
			return false;
		}
	}

	boolean Criteria3_5(TreeMap<String, Integer> controlMap) {
		if (controlMap.size() == 2) {
			return true;
		} else {
			return false;
		}
	}

	boolean Criteria3_6(TreeMap<String, Integer> controlMap) {
		if (numHomozygous(controlMap) > 0) {
			TreeSet<String> controlSet = NewIt.newTreeSet();
			for (String g:controlMap.keySet()) {
				controlSet.add(g.substring(0, 1));
				controlSet.add(g.substring(1, 2));
			}
			if (controlSet.size() == 3) {
				return true;
			} else {
				return false;
			}
		} else {
			return false;
		}
	}

	boolean Criteria3_7(TreeMap<String, Integer> controlMap) {
		if (controlMap.size() == 2) {
			return true;
		} else {
			return false;
		}
	}

	boolean Criteria3_8(TreeMap<String, Integer> controlMap) {
		if (controlMap.size() >= 3) {
			return true;
		} else {
			return false;
		}
	}

	void shuffle(String[] control) {
		// print(control);
		int[] CGSetSize = new int[childrenGenoMap.keySet().size()];
		int index = 0;
		int offset = 0;
		for (String g:childrenGenoMap.keySet()) {
			String geno = g;
			CGSetSize[index] = childrenGenoMap.get(geno).intValue();
			for (int j = 0; j < CGSetSize[index]; j++) {
				control[j + offset] = geno;
			}
			offset += CGSetSize[index];
			index++;
		}
		// print(control);
		int N = control.length;
		for (int i = 0; i < N; i++) {
			int Ind = i + (int) (rnd.nextFloat() * (N - i));
			String tmp = control[i];
			control[i] = control[Ind];
			control[Ind] = tmp;
		}
		// print(control);
	}

	String CompatibleGenotype() {
		TreeMap<String, Integer> alleleMap = new TreeMap<String, Integer>();

		for (String g:childrenGenoMap.keySet()) {
			for (int i = 0; i < 2; i++) {
				if (alleleMap.containsKey(g.substring(i + 0, i + 1))) {
					Integer c = ((Integer) alleleMap.get(g.substring(i + 0, i + 1)));
					c++;
					alleleMap.put(g.substring(i + 0, i + 1), c);
				} else {
					Integer c = new Integer(1);
					alleleMap.put(new String(g.substring(i + 0, i + 1)), c);
				}
			}
		}
		StringBuilder geno = new StringBuilder();

		for (String allele:alleleMap.keySet()) {
			if ((alleleMap.get(allele)).intValue() == 1) {
				geno.append(allele);
			}
		}
		return geno.toString();
	}
}