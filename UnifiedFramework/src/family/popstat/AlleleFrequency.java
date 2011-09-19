package family.popstat;

import statistics.FisherExactTest.FastFisherExactTest;

import org.apache.commons.math.MathException;
import org.apache.commons.math.stat.inference.ChiSquareTest;
import org.apache.commons.math.stat.inference.ChiSquareTestImpl;

public class AlleleFrequency {
	private GenotypeMatrix G;
	private int numMarker;
	private double[][] allelefreq;
	private double[][] genotypefreq;
	private double[][] hw;
	public AlleleFrequency(GenotypeMatrix g) {
		G = g;
		numMarker = g.numMarker;
		allelefreq = new double[numMarker][3];
		genotypefreq = new double[numMarker][4];
		hw = new double[numMarker][2];
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
			hw[i][0] = ff.HDP();
			hw[i][1] = chiHWE(N, NAB, NA);

		}
	}

	public double chiHWE(int N, int NAB, int NA) {
		int n = N;
		int na = NA;
		int nb = 2 * N - na;

		long[] O = { (na - NAB)/2, NAB, (nb - NAB)/2 };
		double n_d = n * 1.0;
		double na_d = na * 1.0;
		double nb_d = nb * 1.0;
		double[] E = { na_d * na_d / (4 * n), na_d * nb_d * 2 / (4*n), nb_d * nb_d / (4*n) };

		ChiSquareTestImpl chi = new ChiSquareTestImpl();
		double p = 0;
		try {
			p = chi.chiSquareTest(E, O);
		} catch (IllegalArgumentException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (MathException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
		return p;
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
			sb.append(" hw.fisher " + String.format("%.3f", hw[i][0]) + " hw.chi " + String.format("%.3f", hw[i][1]));
			sb.append(System.getProperty("line.separator"));
		}
		return sb.toString();
	}
}
