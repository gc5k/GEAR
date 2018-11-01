package gear.util.pop;

import org.apache.commons.math.MathException;
import org.apache.commons.math.random.RandomDataImpl;

import gear.ConstValues;
import gear.family.GenoMatrix.GenotypeMatrix;
import gear.util.Logger;

public class PopStat {
	public static double[][] calAlleleFrequency(GenotypeMatrix G) {
		// [][0]a1 freq; [][1]a2 freq; [][2] missing rate
		// it calculates second allele frequency (so, likely the major one)
		double[][] allelefreq = new double[G.getNumMarker()][3];
		for (int i = 0; i < G.getGRow(); i++) {
			for (int j = 0; j < G.getNumMarker(); j++) {
				int[] c = G.getBiAlleleGenotype(i, j);
				allelefreq[j][c[0]]++;
				allelefreq[j][c[1]]++;
			}
		}

		for (int i = 0; i < G.getNumMarker(); i++) {
			double wa = allelefreq[i][0] + allelefreq[i][1];
			double a = allelefreq[i][0] + allelefreq[i][1] + allelefreq[i][2];
			if (wa > 0) {
				for (int j = 0; j < allelefreq[i].length - 1; j++) {
					allelefreq[i][j] /= wa;
				}
				allelefreq[i][2] /= a;
			} else {
				allelefreq[i][0] = Double.NaN;
				allelefreq[i][0] = Double.NaN;
				allelefreq[i][2] = 1;
			}
		}

		return allelefreq;
	}

	public static double[] calGenoVariance(GenotypeMatrix G) {
		// [][0]a1 freq; [][1]a2 freq; [][2] missing rate
		// it calculates second allele frequency (so, likely the major one)
		double[] axsq = new double[G.getNumMarker()];

		for (int i = 0; i < G.getNumMarker(); i++) {
			int cnt = 0;
			double sq = 0;
			double sm = 0;
			for (int j = 0; j < G.getGRow(); j++) {
				int g = G.getAdditiveScore(j, i);
				if (g != ConstValues.MISSING_GENOTYPE) {
					sq += g * g;
					sm += g;
					cnt++;
				}
			}

			if (cnt > 2) {
				axsq[i] = (sq - cnt * (sm / cnt) * (sm / cnt)) / (cnt - 1);
			}
		}

		return axsq;
	}

	public static double[][] calGenoFrequency(GenotypeMatrix G, boolean isFreq) {
		double[][] gfreq = new double[G.getNumMarker()][3];
		for (int i = 0; i < G.getGRow(); i++) {
			for (int j = 0; j < G.getNumMarker(); j++) {
				int g = G.getAdditiveScoreOnFirstAllele(i, j);
				if (g != ConstValues.MISSING_GENOTYPE) {
					gfreq[j][g]++;
				}
			}
		}

		if (isFreq) {
			for (int i = 0; i < G.getNumMarker(); i++) {
				double wa = gfreq[i][0] + gfreq[i][1] + gfreq[i][2];
				if (wa > 0) {
					gfreq[i][0] /= wa;
					gfreq[i][1] /= wa;
					gfreq[i][2] /= wa;
				}
			}
		}
		return gfreq;
	}

	public static void Imputation(GenotypeMatrix G, boolean isInbred, long seed) {
		Logger.printUserLog("Implementing naive imputatation......");
		double[][] f = calAlleleFrequency(G);
		Logger.printUserLog(
				"Missing genotypes will be imputed according to Hardy-Weinberg proportion for each locus with its estimated allele frequency.");

		RandomDataImpl rnd = new RandomDataImpl();
		rnd.reSeed(seed);

		int cn = 0;
		for (int i = 0; i < G.getNumMarker(); i++) {

			for (int j = 0; j < G.getNumIndivdial(); j++) {
				int genoValue = G.getAdditiveScore(j, i);
				if (genoValue == ConstValues.MISSING_GENOTYPE) {
					cn++;
					int v;
					try {
						if(isInbred) {
							v = rnd.nextBinomial(1, 1 - f[i][0]) * 2;
						} else {
							v = rnd.nextBinomial(2, 1 - f[i][0]);														
						}
						G.setAdditiveScore(j, i, v);
					} catch (MathException e) {
						Logger.handleException(e, "Error in imputation.");
					}
				}
			}
		}
		Logger.printUserLog("In total " + cn + " genotypes have been imputed.");
	}

	public static double[] CalcLDfromDPrime(double[] freq, double[] dprime) {
		double[] D = new double[dprime.length];

		for (int i = 0; i < D.length; i++) {
			if (dprime[i] > 0) {
				D[i] = dprime[i] * Math.min(freq[i] * (1 - freq[i + 1]), freq[i + 1] * (1 - freq[i]));
			} else {
				D[i] = dprime[i] * Math.min(freq[i] * freq[i + 1], (1 - freq[i]) * (1 - freq[i + 1]));
			}
		}

		return D;

	}

	public static double Fst(double F, int n1, int n2, double f1, double f2) {
		double N = n1 + n2;
		double fst = 2 * (n1 / N * (f1 - F) * (f1 - F) + n2 / N * (f2 - F) * (f2 - F))
		/ (F * (1.0D - F));
		return fst;
	}
}
