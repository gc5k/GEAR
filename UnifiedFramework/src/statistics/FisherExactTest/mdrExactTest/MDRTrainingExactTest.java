package statistics.FisherExactTest.mdrExactTest;

import java.util.HashMap;
import java.util.Map.Entry;

import admixture.parameter.Parameter;

import family.mdr.result.Cell;

public class MDRTrainingExactTest {

	private HashMap<String, Cell> classes;
	private double T;
	private int H; // number of high risk groups
	private int L; // number of low risk groups
	private int HPos; // positive individuals in high groups
	private int HNeg; // negative individuals in high groups;
	private int LPos; // positive individuals in low groups;
	private int LNeg; // negative individuals in low groups;
	private int[][] HSub;
	private int[][] LSub;

	private double pbase;
	private double ptruncated;
	private double pobs;
	private double pInt;
	private double pOneTail;

	private int[][] confusion = { { 0, 0 }, { 0, 0 } };

	public MDRTrainingExactTest(HashMap<String, Cell> classes) {
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
//				confusion[0][1] += cell.getPositiveSubjects();
//				confusion[1][0] += cell.getNegativeSubjects();
			}
		}
		HSub = new int[H][2];
		LSub = new int[L][2];

		int c1 = 0;
		int c2 = 0;
		for (Entry<String, Cell> entry : classes.entrySet()) {
			Cell cell = entry.getValue();
			if (cell.getStatus() == 1) {
				HSub[c1][0] = (int) cell.getPositiveSubjects();
				HSub[c1][1] = (int) cell.getNegativeSubjects();
				c1++;
			} else if (cell.getStatus() == 0) {
				LSub[c2][0] = (int) cell.getPositiveSubjects();
				LSub[c2][1] = (int) cell.getNegativeSubjects();
				c2++;
			} else {

			}
		}
		T = 1.0 * (confusion[0][0] + confusion[0][1]) / (confusion[1][0] + confusion[1][1]);
		if (H >= L) {
			for (int i = 0; i < HSub.length; i++) {
				HPos += getNumSubject(true, HSub[i][0], HSub[i][1], T / (1 + T));
			}
			HNeg = confusion[1][0] + confusion[0][0] - HPos;
			LPos = confusion[0][1] + confusion[0][0] - HPos;
			LNeg = confusion[1][1] + confusion[1][0] - HNeg;
		} else {
			for (int i = 0; i < LSub.length; i++) {
				LNeg += getNumSubject(false, LSub[i][0], LSub[i][1], 1.0 / (1 + T));
			}
			LPos = confusion[0][1] + confusion[1][1] - LNeg;
			HNeg = confusion[1][0] + confusion[1][1] - LNeg;
			HPos = confusion[0][0] + confusion[1][0] - HNeg;
		}
	}

	private int getNumSubject(boolean isHigh, int PosSubs, int NegSubs, double T) {
		int n1;
		if (isHigh) {
			n1 = (int) Math.ceil((PosSubs + NegSubs) * T);
			int n2 = PosSubs + NegSubs - n1;

			if (n1 * 1.0 == T * n2 && Parameter.tie == 0) {
				n1++;
				n2--;
			}
		} else {
			n1 = (int) Math.ceil((PosSubs + NegSubs) * T);
			int n2 = PosSubs + NegSubs - n1;

			if (n1 * 1.0 == T * n2 && Parameter.tie == 1) {
				n1--;
				n2++;
			}
		}
		return n1;
	}

	private void base() {

		double p = 0;

		for (int i = 1; i <= HPos + HNeg; i++) {
			p += Math.log(i);
		}
		for (int i = 1; i <= HPos + LPos; i++) {
			p += Math.log(i);
		}
		for (int i = 1; i <= LNeg + HNeg; i++) {
			p += Math.log(i);
		}
		for (int i = 1; i <= LNeg + LPos; i++) {
			p += Math.log(i);
		}

		for (int i = 1; i <= HPos; i++) {
			p -= Math.log(i);
		}
		for (int i = 1; i <= HNeg; i++) {
			p -= Math.log(i);
		}
		for (int i = 1; i <= LPos; i++) {
			p -= Math.log(i);
		}
		for (int i = 1; i <= LNeg; i++) {
			p -= Math.log(i);
		}
		for (int i = 1; i <= HPos + HNeg + LPos + LNeg; i++) {
			p -= Math.log(i);
		}

		pbase = Math.exp(p);
	}

	public double getOneTailP() {
		return pOneTail;
	}

	public void ExactTest() {

		base();

		int upper = LPos < HNeg ? LPos : HNeg;
		upper += HPos - 1;

		int i = HPos;
		int j = 0;

		double p = pbase;

		do {
			ptruncated += p;
			if (i < confusion[0][0]) {
				pInt += p;

			}
			if (i == confusion[0][0]) {
				pobs = p;
			}
			int a = i;
			int b = confusion[0][0] + confusion[0][1] - a;
			int c = confusion[0][0] + confusion[1][0] - a;
			int d = confusion[1][1] + confusion[0][1] - b;
			p = p * b * c / ((a + 1) * (d + 1));
			i++;
			j++;

		} while (i <= upper);

		pOneTail = 1 - pInt / ptruncated;
	}
}
