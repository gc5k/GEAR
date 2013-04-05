package family.RabinowitzLairdAlgorithm.rabinowitz;

import java.util.TreeMap;
import java.util.TreeSet;

import family.RabinowitzLairdAlgorithm.AbstractGenoDistribution;
import gear.util.Logger;
import gear.util.NewIt;

/**
*
* @author Guo-Bo Chen, chenguobo@gmail.com
*/
public class Rabinowitz2 extends AbstractGenoDistribution {

	TreeMap<String, Integer> parentGenoMap;
	String parentgeno1;

	public Rabinowitz2(TreeMap<String, Integer> children, TreeMap<String, Integer> parent) {
		super(children);
		parentGenoMap = new TreeMap<String, Integer>(parent);
		this.parentgeno1 = (String) parentGenoMap.firstKey();
		countChildrenAllele(childrenGenoMap);
		countAllele(childrenGenoMap);
		countParentAllele(parentGenoMap);
		countAllele(parentGenoMap);
	}

	protected void genotypeParents() {
	}

	public String[] getNontransmitted(final String transmitted) {
		return null;
	}

	public String[] getNontransmitted() {
		String[] control = new String[getChildrenNum()];
		TreeMap<String, Integer> controlMap = NewIt.newTreeMap();
		// how to deal with missing data???

		if (childrenGenoMap.size() == 1) {
			if (numHomozygous(childrenGenoMap) == 1 || childrenGenoMap.containsKey(parentgeno1)) {// situation
				// 1
				// System.err.println("Rabinowitz Table2 s1");

				String[] genopool = new String[1];
				genopool[0] = (String) childrenGenoMap.firstKey();
				double[] freq = new double[1];
				freq[0] = 1.0;
				Produce(control, controlMap, genopool, freq);
			} else {// situation 4-1
				// System.err.println("Rabinowitz Table2 s4-1");

				String DeducedGeno1 = ExtractUniqueAllele2Genotype(childrenGenoMap.firstKey(), parentgeno1);
				String[] genopool = new String[2];
				genopool[0] = childrenGenoMap.firstKey();
				genopool[1] = DeducedGeno1;
				double[] freq = new double[2];
				freq[0] = 0.5;
				freq[1] = 1.0;
				Produce(control, controlMap, genopool, freq);
			}
		} else if (childrenGenoMap.size() == 2) {
			if (numHomozygous(childrenGenoMap) == 0) {// situation

				if (childrenGenoMap.containsKey(parentgeno1)) {// 5-1
					// System.err.println("Rabinowitz Table2 s5-1");

					String DeducedGeno1 = ExtractUniqueAllele2Genotype(childrenGenoMap.firstKey(), childrenGenoMap.lastKey());
					String[] genopool = new String[3];
					genopool[0] = childrenGenoMap.firstKey();
					genopool[1] = childrenGenoMap.lastKey();
					genopool[2] = DeducedGeno1;
					double[] freq = new double[3];
					freq[0] = 0.33333;
					freq[1] = 0.66667;
					freq[2] = 1.0;
					do {
						controlMap.clear();
						Produce(control, controlMap, genopool, freq);
					} while (!Criteria2_5(controlMap));
				} else {
					if (numChildrenAllele() == 3) {
						TreeSet<String> CASet = getChildrenAlleleSet();
						if (CASet.contains(parentgeno1.substring(0, 1)) && CASet.contains(parentgeno1.substring(1, 2))) {// 4-2
							// System.err.println("Rabinowitz Table2 s4-2");

							String[] genopool = new String[2];
							genopool[0] = childrenGenoMap.firstKey();
							genopool[1] = childrenGenoMap.lastKey();
							double[] freq = new double[2];
							freq[0] = 0.5;
							freq[1] = 1.0;
							Produce(control, controlMap, genopool, freq);
						} else {// 7-2
							// System.err.println("Rabinowitz Table2 s7-2");

							String[] genopool = new String[4];
							String Deducedgeno1 = ExtractUniqueAllele2Genotype(childrenGenoMap.firstKey(), parentgeno1);
							String Deducedgeno2 = ExtractUniqueAllele2Genotype(childrenGenoMap.lastKey(), parentgeno1);
							genopool[0] = childrenGenoMap.firstKey();
							genopool[1] = childrenGenoMap.lastKey();
							genopool[2] = Deducedgeno1;
							genopool[3] = Deducedgeno2;
							double[] freq = new double[4];
							freq[0] = 0.25;
							freq[1] = 0.5;
							freq[2] = 0.75;
							freq[3] = 1;
							do {
								controlMap.clear();
								Produce(control, controlMap, genopool, freq);
							} while (!Criteria2_7(controlMap));
						}
					} else if (numChildrenAllele() == 4) {// 7-1
						// System.err.println("Rabinowitz Table2 s7-1");

						String[] genopool = new String[4];
						String Deducedgeno1 = ExtractUniqueAllele2Genotype(childrenGenoMap.firstKey(), parentgeno1);
						String Deducedgeno2 = ExtractUniqueAllele2Genotype(childrenGenoMap.lastKey(), parentgeno1);
						genopool[0] = childrenGenoMap.firstKey();
						genopool[1] = childrenGenoMap.lastKey();
						genopool[2] = Deducedgeno1;
						genopool[3] = Deducedgeno2;
						double[] freq = new double[4];
						freq[0] = 0.25;
						freq[1] = 0.5;
						freq[2] = 0.75;
						freq[3] = 1;
						do {
							controlMap.clear();
							Produce(control, controlMap, genopool, freq);
						} while (!Criteria2_7(controlMap));
					}
				}
			}
			if (numHomozygous(childrenGenoMap) == 1) {
				if (childrenGenoMap.containsKey(parentgeno1)) {// 2
					// System.err.println("Rabinowitz Table2 s2");

					shuffle(control);
				} else {// 6-1, 6-2
					// System.err.println("Rabinowitz Table2 s6-1,6-2");

					String[] genopool = new String[4];
					String Deducedgeno1 = ExtractUniqueAllele2Genotype(childrenGenoMap.firstKey(), parentgeno1);
					String Deducedgeno2 = ExtractUniqueAllele2Genotype(childrenGenoMap.lastKey(), parentgeno1);
					genopool[0] = childrenGenoMap.firstKey();
					genopool[1] = childrenGenoMap.lastKey();
					genopool[2] = Deducedgeno1;
					genopool[3] = Deducedgeno2;
					double[] freq = new double[4];
					freq[0] = 0.25;
					freq[1] = 0.5;
					freq[2] = 0.75;
					freq[3] = 1;
					do {
						controlMap.clear();
						Produce(control, controlMap, genopool, freq);
					} while (!Criteria2_6(controlMap));
				}
			}
			if (numHomozygous(childrenGenoMap) == 2) {// 3-1
				// System.err.println("Rabinowitz Table2 s3-1");

				String[] genopool = new String[3];
				String Deducedgeno1 = ExtractUniqueAllele2Genotype(childrenGenoMap.firstKey(), parentgeno1);
				genopool[0] = childrenGenoMap.firstKey();
				genopool[1] = childrenGenoMap.lastKey();
				genopool[2] = Deducedgeno1;
				double[] freq = new double[3];
				freq[0] = 0.25;
				freq[1] = 0.5;
				freq[2] = 1;
				do {
					controlMap.clear();
					Produce(control, controlMap, genopool, freq);
				} while (!Criteria2_3(controlMap));
			}
		} else if (childrenGenoMap.size() == 3) {
			if (numHomozygous(childrenGenoMap) == 0) {
				if (childrenGenoMap.containsKey(parentgeno1)) {// 5-2
					// System.err.println("Rabinowitz Table2 s5-2");

					String[] genopool = new String[3];
					int index = 0;
					for (String g : childrenGenoMap.keySet()) {
						genopool[index] = g;
					}
					double[] freq = new double[3];
					freq[0] = 0.33333;
					freq[1] = 0.66667;
					freq[2] = 1.0;
					do {
						controlMap.clear();
						Produce(control, controlMap, genopool, freq);
					} while (!Criteria2_5(controlMap));
				} else {// 7-3
					// System.err.println("Rabinowitz Table2 s7-3");

					String[] genopool = new String[4];
					int index = 0;
					for (String g : childrenGenoMap.keySet()) {
						genopool[index] = g;
					}
					genopool[3] = CompatibleGenotype();
					double[] freq = new double[4];
					freq[0] = 0.25;
					freq[1] = 0.5;
					freq[2] = 0.75;
					freq[3] = 1;
					do {
						controlMap.clear();
						Produce(control, controlMap, genopool, freq);
					} while (!Criteria2_7(controlMap));
				}
			}
			if (numHomozygous(childrenGenoMap) == 1) {
				if (childrenGenoMap.containsKey(parentgeno1)) {// 6-3, 6-4
					// System.err.println("Rabinowitz Table2 s6-3,6-4");

					String[] genopool = new String[4];
					int index = 0;
					String Deducedgeno1 = new String();
					for (String g : childrenGenoMap.keySet()) {
						genopool[index] = g;
						if (isHeterozygous(genopool[index]) && genopool[index].compareTo(parentgeno1) != 0) {
							Deducedgeno1 = ExtractUniqueAllele2Genotype(parentgeno1, genopool[index]);
						}
						index++;
					}
					genopool[3] = Deducedgeno1;
					double[] freq = new double[4];
					freq[0] = 0.25;
					freq[1] = 0.5;
					freq[2] = 0.75;
					freq[3] = 1;
					do {
						controlMap.clear();
						Produce(control, controlMap, genopool, freq);
					} while (!Criteria2_6(controlMap));
				} else {// 6-5
					// System.err.println("Rabinowitz Table2 s6-5");

					String[] genopool = new String[4];
					int index = 0;
					for (String g : childrenGenoMap.keySet()) {
						genopool[index] = g;
						index++;
					}
					genopool[3] = parentgeno1;
					double[] freq = new double[4];
					freq[0] = 0.25;
					freq[1] = 0.5;
					freq[2] = 0.75;
					freq[3] = 1;
					do {
						controlMap.clear();
						Produce(control, controlMap, genopool, freq);
					} while (!Criteria2_6(controlMap));
				}
			}
			if (numHomozygous(childrenGenoMap) == 2) {// 3-2
				// System.err.println("Rabinowitz Table2 s3-2");

				String[] genopool = new String[3];
				int index = 0;
				for (String g : childrenGenoMap.keySet()) {
					genopool[index] = g;
					index++;
				}
				double[] freq = new double[3];
				freq[0] = 0.25;
				freq[1] = 0.5;
				freq[2] = 1;
				do {
					controlMap.clear();
					Produce(control, controlMap, genopool, freq);
				} while (!Criteria2_3(controlMap));
			}
		} else if (childrenGenoMap.size() == 4) {// 7-4
			// System.err.println("Rabinowitz Table2 s7-4");

			String[] genopool = new String[4];
			int index = 0;
			for (String g : childrenGenoMap.keySet()) {
				genopool[index] = g;
				index++;
			}
			double[] freq = new double[4];
			freq[0] = 0.25;
			freq[1] = 0.5;
			freq[2] = 0.75;
			freq[3] = 1.0;
			do {
				controlMap.clear();
				Produce(control, controlMap, genopool, freq);
			} while (!Criteria2_7(controlMap));
		} else {
			Logger.printUserError("Wrecked in Rabinowitz table 2");
		}
		return control;
	}

	boolean Criteria2_3(TreeMap<String, Integer> controlMap) {
		if (numHomozygous(controlMap) == 2) {
			return true;
		} else {
			return false;
		}
	}

	boolean Criteria2_5(TreeMap<String, Integer> controlMap) {
		if (controlMap.size() > 1 && controlMap.containsKey(parentgeno1)) {
			return true;
		} else {
			return false;
		}
	}

	boolean Criteria2_6(TreeMap<String, Integer> controlMap) {
		if (numHomozygous(controlMap) > 0) {
			for (String g:controlMap.keySet()) {
				if (isHeterozygous(g) && g.compareTo(parentgeno1) != 0) {
					return true;
				}
			}
			return false;
		} else {
			return false;
		}
	}

	boolean Criteria2_7(TreeMap<String, Integer> controlMap) {
		TreeSet<String> controlSet = NewIt.newTreeSet();
		for (String g: controlMap.keySet()) {
			controlSet.add(g.substring(0, 1));
			controlSet.add(g.substring(1, 2));
		}
		if (controlSet.size() == 4) {
			return true;
		} else {
			return false;
		}
	}

	void shuffle(String[] control) {
		int[] CGSetSize = new int[childrenGenoMap.keySet().size()];
		int index = 0;
		int offset = 0;
		for (String geno:childrenGenoMap.keySet()) {
			CGSetSize[index] = ((Integer) childrenGenoMap.get(geno)).intValue();
			for (int j = 0; j < CGSetSize[index]; j++) {
				control[j + offset] = geno;
			}
			offset += CGSetSize[index];
			index++;
		}
		int N = control.length;
		for (int i = 0; i < N; i++) {
			int Ind = i + (int) (rnd.nextFloat() * (N - i));
			String tmp = control[i];
			control[i] = control[Ind];
			control[Ind] = tmp;
		}
	}

	String CompatibleGenotype() {
		TreeMap<String, Integer> alleleMap = NewIt.newTreeMap();

		for (String g:childrenGenoMap.keySet()) {
			for (int i = 0; i < 2; i++) {
				if (alleleMap.containsKey(g.substring(i + 0, i + 1))) {
					Integer c = alleleMap.get(g.substring(i + 0, i + 1));
					c++;
					alleleMap.put(g.substring(i + 0, i + 1), c);
				} else {
					alleleMap.put(new String(g.substring(i + 0, i + 1)), new Integer(1));
				}
			}
		}
		String geno = new String();
		for (String allele:alleleMap.keySet()) {
			if (((Integer) alleleMap.get(allele) == 1)) {
				geno += allele;
			}
		}
		return geno;
	}
}