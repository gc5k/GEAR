package family.popstat;

import statistics.FisherExactTest.FastFisherExactTest;

public class AlleleFrequency {
	private GenotypeMatrix G;
	private int numMarker;
	private double[][] allelefreq;
	private double[][] genotypefreq;
	private double[] hw;
	public AlleleFrequency(GenotypeMatrix g) {
		G = g;
		numMarker = g.numMarker;
		allelefreq = new double[numMarker][3];
		genotypefreq = new double[numMarker][4];
		hw = new double[numMarker];
	}

	public void CalculateAlleleFrequency() {
		int[][] g = G.getG();
		for(int i = 0; i < g.length; i++) {
			for(int j = 0; j < numMarker; j++) {
				int[] c = G.getBiAlleleGenotype(i, j);
				allelefreq[j][c[0]]++; allelefreq[j][c[1]]++;
				int idx = c[0] == 2 ? 3 : (c[0] + c[1]);
				genotypefreq[j][idx]++;
			}
		}

		for(int i = 0; i < numMarker; i++) {
			int N = (int) (genotypefreq[i][0] + genotypefreq[i][1] + genotypefreq[i][2]);
			int NAB = (int) (genotypefreq[i][1]);
			int NA = (int) (allelefreq[i][0] < allelefreq[i][1] ? allelefreq[i][0] : allelefreq[i][1]);
			
			double a = allelefreq[i][0] + allelefreq[i][1] + allelefreq[i][2];
			if (a > 0) {
				for(int j = 0; j < allelefreq[i].length; j++) {
					allelefreq[i][j] /= a;
				}
			} else {
				allelefreq[i][2] = 1;
			}
			double b = genotypefreq[i][0] + genotypefreq[i][1] + genotypefreq[i][2] + genotypefreq[i][3];
			if (b > 0) {
				for(int j = 0; j < genotypefreq[i].length; j++) {
					if (b >= 0) {
						genotypefreq[i][j] /= b;
					}
				}
			} else {
				genotypefreq[i][3] = 1;
			}

			FastFisherExactTest ff = new FastFisherExactTest(N, NAB, NA);
			hw[i] = ff.HDP();
		}
	}
	
	public double[][] getAlleleFrequency() {
		return allelefreq;
	}

	public String toString() {
		StringBuffer sb = new StringBuffer();
		sb.append(G.getG().length);
		sb.append(System.getProperty("line.separator"));
		for(int i = 0; i < allelefreq.length; i++) {
			
			sb.append(String.format("%.3f", allelefreq[i][0]) + " " + String.format("%.3f", allelefreq[i][1]) + " " + String.format("%.3f", allelefreq[i][2])
					+ "; " + String.format("%.3f", genotypefreq[i][0]) 
				+ " " + String.format("%.3f", genotypefreq[i][1]) + " " + String.format("%.3f", genotypefreq[i][2]) + " " + String.format("%.3f", genotypefreq[i][3]));
			sb.append(" hw " + String.format("%.3f", hw[i]));
			sb.append(System.getProperty("line.separator"));
		}
		return sb.toString();
	}
}
