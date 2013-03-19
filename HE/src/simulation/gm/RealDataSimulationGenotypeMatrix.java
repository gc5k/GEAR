package simulation.gm;

import java.util.ArrayList;

import family.pedigree.PersonIndex;
import parameter.Parameter;
import simulation.qc.rowqc.RealDataSimulationQC;

public class RealDataSimulationGenotypeMatrix {

	protected int[][] genotypeMat;
	protected int lenMat;
	protected int lenMatF;
	protected final int shift = 4;
	protected int numMarker = 0;
	protected ArrayList<PersonIndex> pidx;
	
	public RealDataSimulationGenotypeMatrix(RealDataSimulationQC cb) {
		pidx = cb.getSample();
		initial();
	}

	protected void initial() {
		genotypeMat = new int[pidx.size()][];

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
		if (g == 3) {// 01
			return Parameter.INSTANCE.missingGenotype;
		} else {
			return Integer.toString(g);
		}
	}

	public int[][] getG() {
		return genotypeMat;
	}

	public int getGRow() {
		return genotypeMat.length;
	}
	
	public int getGCol() {
		return genotypeMat[0].length;
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
