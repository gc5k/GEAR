package family.mdr.arsenal;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

import org.apache.commons.lang3.ArrayUtils;

import family.mdr.data.PersonIndex;

import util.NewIt;
import admixture.parameter.Parameter;


public class Partition {

	public static Random rnd = new Random(Parameter.seed);

	public static ArrayList<Integer> CVPartition(int n) {
		ArrayList<Integer> g = NewIt.newArrayList();
		for (int i = 0; i < n; i++) {
			g.add(i % Parameter.cv);
		}
		Collections.shuffle(g, rnd);
		return g;
	}

	public static ArrayList<Integer> TTPartition(int n) {
		ArrayList<Integer> g = NewIt.newArrayList();
		double c = Parameter.trgroup;
		for (int i = 0; i < n; i++) {
			if (rnd.nextDouble() <= c) {
				g.add(0);
			} else {
				g.add(1);
			}
		}
		return g;
	}
	
	public static ArrayList<Integer> SexPartition(ArrayList<PersonIndex> PersonTable) {
		ArrayList<Integer> g = NewIt.newArrayList();
		for (PersonIndex pi:PersonTable) {
			int sex = pi.getPerson().getGender();
			if (Parameter.trsex == sex) {
				g.add(0);
			} else {
				g.add(1);
			}
		}
		return g;
	}
	
	public static ArrayList<Integer> ttFilePartition(ArrayList<PersonIndex> PersonTable) {
		ArrayList<Integer> g = NewIt.newArrayList();
		
		for (PersonIndex pi:PersonTable) {
			String Fid = pi.getFamilyID();
			String Iid = pi.getIndividualID();
			
			int idxF = ArrayUtils.indexOf(Parameter.ttArray[0], Fid);
			int idxI = ArrayUtils.indexOf(Parameter.ttArray[1], Iid);
			
			if (idxF>=0 && idxI >= 0 && idxF==idxI) {
				g.add(0);
			} else {
				g.add(1);
			}
		}
		return g;
	}

//	public static ArrayList<Integer> BorderPartition(String fid, String pid, ArrayList<PersonIndex> PersonTable) {
//		ArrayList<Integer> g = NewIt.newArrayList();
//		boolean flag = false;
//		for (PersonIndex pi:PersonTable) {
//			String Fid = pi.getFamilyID();
//			String Iid = pi.getIndividualID();
//			if (!Parameter.reverseborderFlag) {
//				if (flag) {
//					g.add(1);
//				} else {
//					if (Fid.compareTo(fid) == 0 && Iid.compareTo(pid) == 0) {
//						flag = true;
//					}
//					if (!flag) {
//						g.add(0);
//					} else {
//						g.add(1);
//					}
//				}
//			} else {
//				if (flag) {
//					g.add(0);
//				} else {
//					if (Fid.compareTo(fid) == 0 && Iid.compareTo(pid) == 0) {
//						flag = true;
//					}
//					if (flag) {
//						g.add(0);
//					} else {
//						g.add(1);
//					}
//				}
//			}
//		}
//		if (flag == false) {
//			throw new IllegalArgumentException("could not locate individual with option --border.");
//		}
//		return g;
//	}
}
