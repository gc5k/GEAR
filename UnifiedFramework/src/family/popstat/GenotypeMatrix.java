package family.popstat;

import java.util.ArrayList;

import family.mdr.arsenal.MDRConstant;
import family.mdr.data.PersonIndex;
import family.pedigree.design.hierarchy.ChenInterface;

public class GenotypeMatrix {

	protected int[][] genotypeMat;
	protected final int shift = 4;
	protected int numMarker = 0;

	public GenotypeMatrix(ChenInterface cb) {
		initial(cb);
	}
	
	private void initial(ChenInterface cb) {
		ArrayList<PersonIndex> pidx = cb.getSample();
		genotypeMat = new int[pidx.size()][];
		for(int i = 0; i < pidx.size(); i++) {
			genotypeMat[i] = pidx.get(i).getPerson().getAlleleArray();
		}
		numMarker = pidx.get(0).getPerson().getNumMarkers();
	}

	public int getNumMarker() {
		return numMarker;
	}

	public int getGenotypeScore(int idx, int i) {
		int posByte = i >> shift;
		int posBite = (i - (i >> shift << shift)) << 1;
		int g = (genotypeMat[idx][posByte] >> (posBite)) & 3;
		if (g == 1) {//01
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
		int posBite = (i - (i >> shift << shift)) << 1;
		int g = (genotypeMat[idx][posByte] >> posBite) & 3;
		int[] b = {2, 2};
		if (g != 1) {
			b[0] = (g >> 1) & 1;
			b[1] = g & 1;
		}
		return b;
	}
	
	public String getGenotypeScoreString(int idx, int i) {
		int posByte = i >> shift;
		int posBite = (i - (i >> shift << shift)) << 1;
		int g = (genotypeMat[idx][posByte] >> (posBite)) & 3;
		if (g == 1) {//01
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
		for(int i = 0; i < genotypeMat.length; i++) {
			for(int j = 0; j < numMarker; j++) {
				getBiAlleleGenotype(i, j);
			}
		}
		long t2 = System.currentTimeMillis();
		System.err.println(t2-t1);

	}
}
