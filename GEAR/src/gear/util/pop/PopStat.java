package gear.util.pop;

import java.util.ArrayList;
import java.util.stream.IntStream;

import org.apache.commons.math.MathException;
import org.apache.commons.math.random.RandomDataImpl;

import gear.ConstValues;
import gear.data.Person;
import gear.family.GenoMatrix.GenotypeMatrix;
import gear.util.Logger;
import gear.util.NewIt;

public class PopStat {
	
	public static ArrayList<ArrayList<Integer>> punchMissingGenoMT(GenotypeMatrix G, int threadNum) {
		final ArrayList<ArrayList<Integer>> missList = NewIt.newArrayList();

		for (int i = 0; i < G.getNumIndivdial(); i++) {
			ArrayList<Integer> mL = NewIt.newArrayList();
			missList.add(mL);
		}

		int markCnt = G.getNumMarker();
		final int cpus = threadNum < markCnt ? threadNum : 1;

		Thread[] computeThreads = new Thread[cpus];
		final int[] taskProgresses = new int[cpus];
		final int taskSize = markCnt / cpus;

		for (int c = 0; c < cpus; ++c) {
			final int threadIndex = c;
			Thread thread = new Thread() {
				public void run() {
					int markStart = threadIndex * taskSize;
					int markEnd = threadIndex < (cpus - 1) ? (taskSize * (threadIndex+1)) : markCnt;
					
					int taskProgress = 0;
					for (int i = markStart; i < markEnd; i++) {
						taskProgresses[threadIndex] = ++taskProgress;

						for (int j = 0; j < G.getNumIndivdial(); j++) {
							int g = G.getAdditiveScore(j,i);
							if (g == Person.MissingGenotypeCode) {
								ArrayList<Integer> mL = missList.get(j);
								mL.add(i);
								missList.set(j, mL);
							}
						}
					}
				}
			};
			thread.start();
			computeThreads[c] = thread;
		}
		
		Thread progressDisplayThread = new Thread() {
			public void run() {
				int totalProgress;
				do {
					totalProgress = IntStream.of(taskProgresses).sum();
					float percentage = Math.min(100f, (float) totalProgress / G.getNumMarker() * 100f);
					System.out.print(String.format(
							"\r[INFO] Calculating missing genotype punchcard, %.2f%% completed...",
							percentage));
					try {
						Thread.sleep(1000);
					} catch (InterruptedException e) {
					}
				} while (totalProgress < G.getNumMarker());
				System.out.println("Done.");
				System.out.print("");
			}
		};
		progressDisplayThread.start();

		for (int i = 0; i < computeThreads.length; ++i) {
			try {
				computeThreads[i].join();
			} catch (InterruptedException e) {
				Logger.handleException(e, String.format("Compute thread %d is interrupted.", i));
			}
		}
		try {
			progressDisplayThread.join();
		} catch (InterruptedException e) {
			Logger.printUserError("Progress display thread is interrupted.");
		}

		return missList;
	}


	public static ArrayList<ArrayList<Integer>> punchMissingGeno(GenotypeMatrix G) {
		ArrayList<ArrayList<Integer>> missList = NewIt.newArrayList();

		for (int i = 0; i < G.getNumIndivdial(); i++) {
			ArrayList<Integer> mL = NewIt.newArrayList();
			missList.add(mL);
		}

		for (int i = 0; i < G.getNumMarker(); i++) {
			for (int j = 0; j < G.getNumIndivdial(); j++) {
				int g = G.getAdditiveScore(j,i);
				if (g == Person.MissingGenotypeCode) {
					ArrayList<Integer> mL = missList.get(j);
					mL.add(i);
//					missList.set(j, mL);
				}
			}
		}
		
		return missList;
	}

	
	public static float[][] calLocusStatMT(GenotypeMatrix G, int threadNum) {
		final float[][] axsq = new float[G.getNumMarker()][4];

		int markCnt = G.getNumMarker();
		final int cpus = threadNum < markCnt ? threadNum : 1;

		Thread[] computeThreads = new Thread[cpus];
		final int[] taskProgresses = new int[cpus];
		final int taskSize = markCnt / cpus;

		for (int c = 0; c < cpus; ++c) {
			final int threadIndex = c;
			Thread thread = new Thread() {
				public void run() {

					int markStart = threadIndex * taskSize;
					int markEnd = threadIndex < (cpus - 1) ? (taskSize * (threadIndex+1)) : markCnt;
					
					int taskProgress = 0;
					for (int i = markStart; i < markEnd; i++) {

						taskProgresses[threadIndex] = ++taskProgress;

						float cnt = 0;
						float sq = 0;
						float sm = 0;
						for (int j = 0; j < G.getNumIndivdial(); j++) {
							int g = G.getAdditiveScore(j, i);
							if (g != ConstValues.MISSING_GENOTYPE) {
								sq += g * g;
								sm += g;
								cnt++;
							}
						}

						if (cnt == 0.0) {
							axsq[i][0] = Float.NaN;
							axsq[i][1] = Float.NaN;
							axsq[i][2] = Float.NaN;
							axsq[i][3] = 1;
							continue;
						}
						axsq[i][0] = 1 - sm / cnt /2;
						axsq[i][1] = sm / cnt /2;

						axsq[i][3] =  (G.getNumIndivdial()-cnt)/(G.getNumIndivdial());
						if (cnt > 2) {
							axsq[i][2] = (sq - cnt * (sm / cnt) * (sm / cnt)) / (cnt - 1);
						}
					}
				}
			};
			thread.start();
			computeThreads[c] = thread;
		}

		Thread progressDisplayThread = new Thread() {
			public void run() {
				int totalProgress;
				do {
					totalProgress = IntStream.of(taskProgresses).sum();
					float percentage = Math.min(100f, (float) totalProgress / G.getNumMarker() * 100f);
					System.out.print(String.format(
							"\r[INFO] Calculating locus stats, %.2f%% completed...",
							percentage));
					try {
						Thread.sleep(1000);
					} catch (InterruptedException e) {
					}
				} while (totalProgress < G.getNumMarker());
				System.out.println("Done.");
				System.out.print("");
			}
		};
		progressDisplayThread.start();

		for (int i = 0; i < computeThreads.length; ++i) {
			try {
				computeThreads[i].join();
			} catch (InterruptedException e) {
				Logger.handleException(e, String.format("Compute thread %d is interrupted.", i));
			}
		}
		try {
			progressDisplayThread.join();
		} catch (InterruptedException e) {
			Logger.printUserError("Progress display thread is interrupted.");
		}
		
		return axsq;
	}

	public static float[][] calLocusStat(GenotypeMatrix G) {
		// [][0]a1 freq; [][1]a2 var; [][2] missing rate
		// it calculates second allele frequency (so, likely the major one)
		final float[][] axsq = new float[G.getNumMarker()][4];

		for (int i = 0; i < G.getNumMarker(); i++) {
			float cnt = 0;
			float sq = 0;
			float sm = 0;
			for (int j = 0; j < G.getNumIndivdial(); j++) {
				int g = G.getAdditiveScore(j, i);
				if (g != ConstValues.MISSING_GENOTYPE) {
					sq += g * g;
					sm += g;
					cnt++;
				}
			}

			if (cnt == 0.0) {
				axsq[i][0] = Float.NaN;
				axsq[i][1] = Float.NaN;
				axsq[i][2] = Float.NaN;
				axsq[i][3] = 1;
				continue;
			}
			axsq[i][0] = 1 - sm / cnt /2;
			axsq[i][1] = sm / cnt /2;

			axsq[i][3] =  (G.getNumIndivdial()-cnt)/(G.getNumIndivdial());
			if (cnt > 2) {
				axsq[i][2] = (sq - cnt * (sm / cnt) * (sm / cnt)) / (cnt - 1);
			}
		}

		return axsq;
	}


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
				int g = G.getAdditiveScore(i, j);
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
