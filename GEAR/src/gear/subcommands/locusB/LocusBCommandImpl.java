package gear.subcommands.locusB;

import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.stream.IntStream;

import gear.family.pedigree.Hukou;
import gear.family.pedigree.file.BEDReader;
import gear.family.pedigree.file.MapFile;
import gear.family.pedigree.file.SNP;
import gear.family.plink.PLINKBinaryParser;
import gear.family.plink.PLINKParser;
import gear.qc.sampleqc.SampleFilter;
import gear.qc.snpqc.SNPFilter;
import gear.qc.snpqc.SNPFilterPostQC;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.FileUtil;
import gear.util.Logger;

public class LocusBCommandImpl extends CommandImpl {
	private LocusBCommandArguments locusArgs;
	private MapFile map;
	private SNPFilter snpFilter;
	private BEDReader bed;
	private PrintStream resultFile;
	private DecimalFormat fmt1 = new DecimalFormat("0.0000");
	private DecimalFormat fmt2 = new DecimalFormat("0.00E000");

	@Override
	public void execute(CommandArguments cmdArgs) {
		
		locusArgs = (LocusBCommandArguments) cmdArgs;

		PLINKParser pp = PLINKParser.create(locusArgs);
		pp.parseSmallFiles();
		map = pp.getMapData();
		snpFilter = pp.getSNPFilter();
		SampleFilter samFilter = new SampleFilter(pp.getPedigreeData(), locusArgs);
		samFilter.qualification();

		resultFile = FileUtil.CreatePrintStream(locusArgs.getOutRoot() + ".locus");
		resultFile.println("SNP\tCHR\tBP\tRefAllele\tAltAllele\tFreq\tVar\tEVar\tAA\tAa\taa\tnChr");

		bed = (BEDReader) pp.getPedigreeData();
		if (!locusArgs.isThreadNum() || (locusArgs.isThreadNum() && locusArgs.getThreadNum() == 1)) {
			executeBedSnpMajor();
		} else {
			executeBedSnpMajorMT();
		}
		return;
	}

	private void executeBedSnpMajor() {
		int numMarkers = map.getMarkerNumberOriginal();
		int numSamples = bed.getNumIndividuals();
		int workingSnpIndex = 0;
		SNPFilterPostQC snpPostQC = new SNPFilterPostQC(locusArgs);
		
		ArrayList<Hukou> hkBook = bed.getHukouBook();
		
		byte[] genotypes = new byte[bed.getByteCountPerRow()];

		for (int i = 0; i < numMarkers; ++i) {
			if (snpFilter.isSnpIncluded(i)) {
				int genoCnt_AA = 0;
				int genoCnt_Aa = 0;
				int genoCnt_aa = 0;
				int missingCnt = 0;
				int sum = 0;
				int squareSum = 0;
				bed.read(genotypes);
				for (int byteIndex = 0, sampleIndex = 0; byteIndex < genotypes.length; ++byteIndex) {
					byte nextByte = genotypes[byteIndex];
					for (int k = 0; k < 8 && sampleIndex < numSamples; k += 2, ++sampleIndex) {
						if (!hkBook.get(sampleIndex).isAvailable()) {
							continue;
						}

						int genotype = (nextByte >> k) & 0b11;

						switch (genotype) {
						case PLINKBinaryParser.HOMOZYGOTE_FIRST:
							++genoCnt_AA;
							break;
						case PLINKBinaryParser.HETEROZYGOTE:
							++genoCnt_Aa;
							sum += 1;
							squareSum += 1;
							break;
						case PLINKBinaryParser.HOMOZYGOTE_SECOND:
							++genoCnt_aa;
							sum += 2;
							squareSum += 4;
							break;
						case PLINKBinaryParser.MISSING_GENOTYPE:
							++missingCnt;
							break;
						}
					}
				}
				int validSampleCnt = numSamples - missingCnt;
				float variance = 0;
				float alleleFreq0 = 0;
				float alleleFreq1 = 0;

				
				if (validSampleCnt > 2) {
					float average = (float)sum / validSampleCnt;
					alleleFreq1 = average/2;
					alleleFreq0 = 1 - alleleFreq1;
					variance = (squareSum - validSampleCnt * average * average) / (validSampleCnt - 1);
				}
				float eVariance = calculateEVariance(alleleFreq0, alleleFreq1);

				boolean isPassPostQC = snpPostQC.isPassPostQC(alleleFreq0, missingCnt * 1.0D / numSamples);
				if (!isPassPostQC) continue;

				SNP snp = map.getSNP(workingSnpIndex++);
				printResultOfSNP(snp, alleleFreq0, variance, eVariance, genoCnt_AA, genoCnt_Aa, genoCnt_aa);
			} else {
				bed.skipOneRow();
			}
		}
		resultFile.close();
		snpPostQC.printPostQCSummary();
		Logger.printUserLog("Saved " + workingSnpIndex + " results to " + locusArgs.getOutRoot() + ".locus.");
	}

	private void executeBedSnpMajorMT() {
		int numMarkers = map.getMarkerNumberOriginal();
		int numSamples = bed.getNumIndividuals();
		int workingSnpIndex = 0;
		SNPFilterPostQC snpPostQC = new SNPFilterPostQC(locusArgs);

		final ArrayList<Hukou> hkBook = bed.getHukouBook();
		final float[][] locusStat = new float[numMarkers][6];

		int cpus = locusArgs.isThreadNum()? locusArgs.getThreadNum():1;
		Logger.printUserLog("Calculating locus statistics with " + cpus + " threads.");
		Thread[] computeThreads = new Thread[cpus];
		final int[] taskProgresses = new int[cpus];
		final int smallTaskSize = numMarkers / cpus;
		final int largeTaskSize = numMarkers - smallTaskSize * (cpus - 1);

		final byte[][] genotypes = new byte[cpus][];

		for (int i = 0; i < cpus; i++) {
			byte[] g;
			if (i < (cpus -1)) {
				g = new byte[bed.getByteCountPerRow()*smallTaskSize];
			} else {
				g = new byte[bed.getByteCountPerRow()*largeTaskSize];
			}
			bed.read(g);
			genotypes[i] = g;
		}
		
		for (int c = 0; c < cpus; ++c) {
			final int threadIndex = c;
			Thread thread = new Thread() {
				public void run() {
					int markerStart = (threadIndex) * smallTaskSize;
					int markerCnt = threadIndex < (cpus-1) ? smallTaskSize:largeTaskSize;
					int taskProgress = 0;
					for (int i = 0; i < markerCnt; i++) {
						taskProgresses[threadIndex] = ++taskProgress;

						if (snpFilter.isSnpIncluded(i+markerStart)) {
							int genoCnt_AA = 0;
							int genoCnt_Aa = 0;
							int genoCnt_aa = 0;
							int missingCnt = 0;
							int sum = 0;
							int squareSum = 0;
							for (int byteIndex = 0, sampleIndex = 0; byteIndex < bed.getByteCountPerRow(); ++byteIndex) {
								byte nextByte = genotypes[threadIndex][byteIndex + i * bed.getByteCountPerRow()];
								for (int k = 0; k < 8 && sampleIndex < numSamples; k += 2, ++sampleIndex) {
									if (!hkBook.get(sampleIndex).isAvailable()) {
										continue;
									}

									int genotype = (nextByte >> k) & 0b11;

									switch (genotype) {
									case PLINKBinaryParser.HOMOZYGOTE_FIRST:
										++genoCnt_AA;
										break;
									case PLINKBinaryParser.HETEROZYGOTE:
										++genoCnt_Aa;
										sum += 1;
										squareSum += 1;
										break;
									case PLINKBinaryParser.HOMOZYGOTE_SECOND:
										++genoCnt_aa;
										sum += 2;
										squareSum += 4;
										break;
									case PLINKBinaryParser.MISSING_GENOTYPE:
										++missingCnt;
										break;
									}
								}
							}
							int validSampleCnt = numSamples - missingCnt;
							float variance = 0;
							float alleleFreq0 = 0;
							float alleleFreq1 = 0;

							if (validSampleCnt > 2) {
								float average = (float)sum / validSampleCnt;
								alleleFreq1 = average/2;
								alleleFreq0 = 1 - alleleFreq1;
								variance = (squareSum - validSampleCnt * average * average) / (validSampleCnt - 1);
							}

							float eVariance = calculateEVariance(alleleFreq0, alleleFreq1);
							locusStat[i+markerStart][0] = alleleFreq0;
							locusStat[i+markerStart][1] = variance;
							locusStat[i+markerStart][2] = eVariance;
							locusStat[i+markerStart][3] = genoCnt_AA;
							locusStat[i+markerStart][4] = genoCnt_Aa;
							locusStat[i+markerStart][5] = genoCnt_aa;
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
					float percentage = Math.min(100f, (float) totalProgress / numMarkers * 100f);
					System.out.print(String.format(
							"\r[INFO] Calculating locus stats, %.2f%% completed...",
							percentage));
					try {
						Thread.sleep(1000);
					} catch (InterruptedException e) {
					}
				} while (totalProgress < numMarkers);
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

		for(int i = 0; i < numMarkers; i++) {
			if (snpFilter.isSnpIncluded(i)) {
				boolean isPassPostQC = snpPostQC.isPassPostQC(locusStat[i][0], locusStat[i][5] * 1.0D / numSamples);
				if (!isPassPostQC) continue;

				SNP snp = map.getSNP(i);
				printResultOfSNP(snp, locusStat[i][0], locusStat[i][1], locusStat[i][2],
						(int) locusStat[i][3], (int) locusStat[i][4], (int) locusStat[i][5]);
				workingSnpIndex++;
			}
		}

		Logger.printUserLog("Saved " + workingSnpIndex + " results to " + locusArgs.getOutRoot() + ".locus.");
	}
	
	private float calculateEVariance(float alleleFreq0, float alleleFreq1) {
		return locusArgs.isInbred() ? 4 * alleleFreq0 * alleleFreq1 : 2 * alleleFreq0 * alleleFreq1;
	}
	
	private void printResultOfSNP(
			SNP snp,
			float alleleFreq0,
			float alleleVar,
			float eVar,
			int genoCnt_AA,
			int genoCnt_Aa,
			int genoCnt_aa) {
		if (Float.isNaN(alleleFreq0)) {
			resultFile.println(snp.getName() + "\t" + snp.getChromosome() + "\t" + snp.getPosition() + "\t"
					+ snp.getFirstAllele() + "\t" + snp.getSecAllele() + "\t"
					+ "NA" + "\t"
					+ "NA" + "\t"
					+ "NA" + "\t"
					+ genoCnt_AA + "\t"
					+ genoCnt_Aa + "\t"
					+ genoCnt_aa + "\t"
					+ ((genoCnt_AA + genoCnt_Aa + genoCnt_aa) << 1));			
		} else {
			resultFile.println(snp.getName() + "\t" + snp.getChromosome() + "\t" + snp.getPosition() + "\t"
					+ snp.getFirstAllele() + "\t" + snp.getSecAllele() + "\t"
					+ (alleleFreq0 > 0.0001 ? fmt1.format(alleleFreq0) : fmt2.format(alleleFreq0)) + "\t"
					+ (alleleVar > 0.0001 ? fmt1.format(alleleVar) : fmt2.format(alleleVar)) + "\t"
					+ (eVar > 0.001 ? fmt1.format(eVar) : fmt2.format(eVar)) + "\t"
					+ genoCnt_AA + "\t"
					+ genoCnt_Aa + "\t"
					+ genoCnt_aa + "\t"
					+ ((genoCnt_AA + genoCnt_Aa + genoCnt_aa) << 1));			
		}
	}
}
