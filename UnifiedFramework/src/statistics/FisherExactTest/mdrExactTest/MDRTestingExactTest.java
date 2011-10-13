package statistics.FisherExactTest.mdrExactTest;

import java.util.HashMap;
import java.util.Map.Entry;
import family.mdr.result.Cell;

public class MDRTestingExactTest {

	private HashMap<String, Cell> classes;

	private int H; // number of high risk groups
	private int L; // number of low risk groups
//	private int[][] HSub;
//	private int[][] LSub;

	private double pbase;
	private double pOneTail;

	private int[][] confusion = { { 0, 0 }, { 0, 0 } };

	// definition of the confusion table
	//       | Hgroup   | Lgroup
	// ---------------------------
	// HGeno | a (HPos) | b (LPos)|
	// --------------------------
	// LGeno | c (HNeg) | d (LNeg)|
	// --------------------------

	public MDRTestingExactTest(HashMap<String, Cell> classes) {
		this.classes = classes;
		initial();
		ExactTest();
	}

	private void initial() {
		for (Entry<String, Cell> entry : classes.entrySet()) {
			Cell cell = entry.getValue();
			if (cell.getStatus() == 1) {
				H++;
				confusion[0][0] += cell.getPositiveSubjects();
				confusion[1][0] += cell.getNegativeSubjects();
			} else if (cell.getStatus() == 0) {
				L++;
				confusion[0][1] += cell.getPositiveSubjects();
				confusion[1][1] += cell.getNegativeSubjects();
			} else {
				confusion[0][1] += cell.getPositiveSubjects();
				confusion[1][0] += cell.getNegativeSubjects();
			}
		}
//		HSub = new int[H][2];
//		LSub = new int[L][2];
//
//		int c1 = 0;
//		int c2 = 0;
//		for (Entry<String, Cell> entry : classes.entrySet()) {
//			Cell cell = entry.getValue();
//			if (cell.getStatus() == 1) {
//				HSub[c1][0] = (int) cell.getPositiveSubjects();
//				HSub[c1][1] = (int) cell.getNegativeSubjects();
//				c1++;
//			} else if (cell.getStatus() == 0) {
//				LSub[c2][0] = (int) cell.getPositiveSubjects();
//				LSub[c2][1] = (int) cell.getNegativeSubjects();
//				c2++;
//			} else {
//				
//			}
//		}
	}

	private void base() {

		double p = 0;

		for (int i = 1; i <= confusion[0][0] + confusion[1][0]; i++) {
			p += Math.log(i);
		}
		for (int i = 1; i <= confusion[0][0] + confusion[0][1]; i++) {
			p += Math.log(i);
		}
		for (int i = 1; i <= confusion[1][0] + confusion[1][1]; i++) {
			p += Math.log(i);
		}
		for (int i = 1; i <= confusion[1][1] + confusion[0][1]; i++) {
			p += Math.log(i);
		}

		for (int i = 1; i <= confusion[0][0]; i++) {
			p -= Math.log(i);
		}
		for (int i = 1; i <= confusion[0][1]; i++) {
			p -= Math.log(i);
		}
		for (int i = 1; i <= confusion[1][0]; i++) {
			p -= Math.log(i);
		}
		for (int i = 1; i <= confusion[1][1]; i++) {
			p -= Math.log(i);
		}
		int N = confusion[0][0] + confusion[0][1] + confusion[1][0] + confusion[1][1];
		for (int i = 1; i <= N; i++) {
			p -= Math.log(i);
		}
		pbase = Math.exp(p);
	}

	public double getOneTailP() {
		return pOneTail;
	}

	public void ExactTest() {

		base();

		int i = confusion[0][0];
		double p = pbase;

		if (confusion[0][0] > confusion[0][1]) {
			int upper = confusion[0][1] < confusion[1][0] ? confusion[0][1] : confusion[1][0];
			upper += confusion[0][0] - 1;

			do {
				pOneTail += p;
				int a = i;
				int b = confusion[0][0] + confusion[0][1] - a;
				int c = confusion[0][0] + confusion[1][0] - a;
				int d = confusion[1][1] + confusion[0][1] - b;
				p = p * b * c / ((a + 1) * (d + 1));
				i++;
			} while (i <= upper);
		} else {
			int lower = 0 > (confusion[0][0] - confusion[1][1]) ? 0 : confusion[0][0] - confusion[1][1];

			do {
				pOneTail += p;
				int a = i;
				int b = confusion[0][0] + confusion[0][1] - a;
				int c = confusion[0][0] + confusion[1][0] - a;
				int d = confusion[1][1] + confusion[0][1] - b;
				p = p * (a * d) / ((b + 1) * (c + 1)) ;
				i--;
			} while (i >= lower);
			pOneTail -= pbase;
			pOneTail = 1 - pOneTail;

		}
	}
}
