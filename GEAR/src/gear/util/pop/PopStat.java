package gear.util.pop;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.IntStream;

import org.apache.commons.math.MathException;
import org.apache.commons.math.random.RandomDataImpl;

import gear.ConstValues;
import gear.data.Person;
import gear.family.GenoMatrix.GenotypeMatrix;
import gear.util.Logger;
import gear.util.NewIt;

public class PopStat {
	
	public static ArrayList<List<Integer>> punchMissingGenoMT(GenotypeMatrix G, int threadNum) {
		final ArrayList<List<Integer>> missList = NewIt.newArrayList();

		for (int i = 0; i < G.getNumIndivdial(); i++) {
			List<Integer> mL = Collections.synchronizedList(NewIt.newArrayList());
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
								List<Integer> mL = missList.get(j);
								mL.add(i);
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
		
		for (List<Integer> missingMarkersOfOneSample : missList) {
			Collections.sort(missingMarkersOfOneSample);
		}

		return missList;
	}

	public static int[][] punchMissingGenoMT2(GenotypeMatrix G, int threadNum) {
		final int[][] missList = new int[G.getNumIndivdial()][];

		int indCnt = G.getNumIndivdial();
		final int cpus = threadNum < indCnt ? threadNum : 1;

		Logger.printUserLog("Calculating missing genotype punchcard with " + cpus + " threads.");

		Thread[] computeThreads = new Thread[cpus];
		final int[] taskProgresses = new int[cpus];
		final int taskSize = indCnt / cpus;

		for (int c = 0; c < cpus; ++c) {
			final int threadIndex = c;
			Thread thread = new Thread() {
				public void run() {
					int indStart = threadIndex * taskSize;
					int indEnd = threadIndex < (cpus - 1) ? (taskSize * (threadIndex+1)) : indCnt;

					int taskProgress = 0;
					for (int i = indStart; i < indEnd; i++) {
						taskProgresses[threadIndex] = ++taskProgress;
						ArrayList<Integer> mList = NewIt.newArrayList();

						for (int j = 0; j < G.getNumMarker(); j++) {
							int g = G.getAdditiveScore(i,j);
							if (g == Person.MissingGenotypeCode) {
								mList.add(j);
							}
						}
						if (mList.size() > 0) {
							missList[i] = new int[mList.size()];
							for (int k = 0; k < mList.size(); k++) {
								missList[i][k] = mList.get(k).intValue();
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
					float percentage = Math.min(100f, (float) totalProgress / G.getNumIndivdial() * 100f);
					System.out.print(String.format(
							"\r[INFO] Calculating missing genotype punchcard, %.2f%% completed...",
							percentage));
					try {
						Thread.sleep(1000);
					} catch (InterruptedException e) {
					}
				} while (totalProgress < G.getNumIndivdial());
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

	public static ArrayList<List<Integer>> punchMissingGeno(GenotypeMatrix G) {
		ArrayList<List<Integer>> missList = NewIt.newArrayList();

		for (int i = 0; i < G.getNumIndivdial(); i++) {
			ArrayList<Integer> mL = NewIt.newArrayList();
			missList.add(mL);
		}

		for (int i = 0; i < G.getNumMarker(); i++) {
			for (int j = 0; j < G.getNumIndivdial(); j++) {
				int g = G.getAdditiveScore(j,i);
				if (g == Person.MissingGenotypeCode) {
					List<Integer> mL = missList.get(j);
					mL.add(i);
				}
			}
		}
		
		return missList;
	}

	public static float[][] calLocusStatFullMT(GenotypeMatrix G, int threadNum) {
		final float[][] axsq = new float[G.getNumMarker()][7];

		int markCnt = G.getNumMarker();
		final int cpus = threadNum < markCnt ? threadNum : 1;

		Logger.printUserLog("Calculating locus statistics with " + cpus + " threads.");
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
						float[] genoCnt = new float[3];

						for (int j = 0; j < G.getNumIndivdial(); j++) {
							int g = G.getAdditiveScore(j, i);
							if (g != ConstValues.MISSING_GENOTYPE) {
								sq += g * g;
								sm += g;
								cnt++;
								genoCnt[g]++;
							}
						}

						if (cnt == 0.0) {
							axsq[i][0] = Float.NaN;
							axsq[i][1] = Float.NaN;

							axsq[i][2] = Float.NaN;

							axsq[i][3] = Float.NaN;
							axsq[i][4] = Float.NaN;
							axsq[i][5] = Float.NaN;
							axsq[i][6] = 1;
							continue;
						}
						axsq[i][0] = 1 - sm / cnt /2;
						axsq[i][1] = sm / cnt /2;

						if (cnt > 2) {
							axsq[i][2] = (sq - cnt * (sm / cnt) * (sm / cnt)) / (cnt - 1);
						}
						axsq[i][3] = genoCnt[0];
						axsq[i][4] = genoCnt[1];
						axsq[i][5] = genoCnt[2];

						axsq[i][6] =  (G.getNumIndivdial()-cnt)/(G.getNumIndivdial());

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
	
	public static double FstWC(double F, int n1, int n2, double f1, double f2) {
		double nn= 1.0*(n1*n2)/(n1+n2);
		double h1 = f1*(1-f1);
		double h2 = f2*(1-f2);
		double dsq = (f1-f2)*(f1-f2);
		double d1= (2 * nn * (1/(n1 * 1.0d + n2 * 1.0d-2))*(n1 * 1.0d * h1 + n2 * 1.0d * h2));
		double d2= (nn * dsq+(2*nn-1)*1/(n1+n2-2)*(n1*h1+n2*h2));
		double Fwc = 1.0 - d1/d2;
		return Fwc;
	}

	public static double FstNei(double F, int n1, int n2, double f1, double f2) {
		double p = (f1+f2)/2;
		double q = 1 - p;
		double dsq = (f1 - f2)*(f1 - f2);
		double FNei = dsq/(2*p*q);
		return FNei;
	}

	public static double FstHudson(double F, int n1, int n2, double f1, double f2) {
		double h1 = f1*(1-f1);
		double h2 = f2*(1-f2);
		double dsq = (f1 - f2)*(f1 - f2);
		double FHudson = (dsq-h1/(n1-1)-h2/(n2-1))/(h1+h2);
		return FHudson;
	}
}
