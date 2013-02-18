package family.popstat;

import java.util.ArrayList;

import family.mdr.arsenal.MDRConstant;
import family.mdr.data.PersonIndex;
import family.pedigree.design.hierarchy.AJHG2008;
import family.pedigree.design.hierarchy.ChenInterface;

public class GenotypeMatrix {

	protected int[][] genotypeMat;
	protected int lenMat;
	protected int lenMatF;
	protected final int shift = 4;
	protected int numMarker = 0;
	protected ArrayList<PersonIndex> pidx;
	protected boolean isAJHG2008;
	
	public GenotypeMatrix(ChenInterface cb) {
		pidx = cb.getSample();
		isAJHG2008 = cb instanceof AJHG2008;
		initial();
	}

	protected void initial() {
		if (isAJHG2008) {
			genotypeMat = new int[pidx.size()/2][];
		} else {
			genotypeMat = new int[pidx.size()][];
		}

		int c1 = 0;
		for (PersonIndex pi : pidx) {
			if (!pi.isPseudo()) {
				genotypeMat[c1++] = pi.getPerson().getAlleleArray();
			}
		}
		numMarker = pidx.get(0).getPerson().getNumMarkers();

	}

	public int getNumMarker() {
		return numMarker;
	}

	public int getGenotypeScore(int idx, int i) {
		int posByte = i >> shift;
		int posBite = (i & 0xf) << 1;
		int g = (genotypeMat[idx][posByte] >> (posBite)) & 3;
		if (g == 1) {// 01
			return 2;
		} else {
			if (g == 2) {
				return 1;
			} else {
				return g;
			}
		}
	}

	public int[] getBiAlleleGenotype(int idx, int i) {
		int posByte = i >> shift;
		int posBite = (i & 0xf) << 1;
		int g = (genotypeMat[idx][posByte] >> posBite) & 3;
		int[] b = { 2, 2 };
		if (g != 1) {
			b[0] = (g >> 1) & 1;
			b[1] = g & 1;
		}
		return b;
	}

	public String getGenotypeScoreString(int idx, int i) {
		int posByte = i >> shift;
		int posBite = (i & 0xf) << 1;
		int g = (genotypeMat[idx][posByte] >> (posBite)) & 3;
		if (g == 1) {// 01
			return MDRConstant.missingGenotype;
		} else {
			if (g == 2) {
				return Integer.toString(1);
			} else {
				return Integer.toString(g);
			}
		}
	}

	public int[][] getG() {
		return genotypeMat;
	}

	public void Test() {
		long t1 = System.currentTimeMillis();
		System.err.println(t1);
		for (int i = 0; i < genotypeMat.length; i++) {
			for (int j = 0; j < numMarker; j++) {
				getBiAlleleGenotype(i, j);
			}
		}
		long t2 = System.currentTimeMillis();
		System.err.println(t2 - t1);

	}
}
