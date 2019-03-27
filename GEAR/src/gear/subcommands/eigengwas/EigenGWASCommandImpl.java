package gear.subcommands.eigengwas;

import gear.ConstValues;
import gear.data.InputDataSet2;
import gear.family.GenoMatrix.GenotypeMatrix;
import gear.family.pedigree.file.SNP;
import gear.family.plink.PLINKParser;
import gear.qc.sampleqc.SampleFilter;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.stream.IntStream;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ChiSquaredDistributionImpl;
import org.apache.commons.math.stat.regression.SimpleRegression;

public class EigenGWASCommandImpl extends CommandImpl {
	private EigenGWASCommandArguments eigenArgs;
	private SampleFilter sf;
	private GenotypeMatrix pGM;

	private int traitIdx;
	private InputDataSet2 data = new InputDataSet2();
	private EigenGWASResult[] eGWASResult = null;

	private double lambdaGC = 1;
	private double[] GC;
	private double[] pQuantile;

	public void execute(CommandArguments cmdArgs) {

		eigenArgs = (EigenGWASCommandArguments) cmdArgs;

		PLINKParser pp = PLINKParser.parse(cmdArgs);

		traitIdx = eigenArgs.getSelectedPhenotype(0);
		data.addFile(eigenArgs.getFam()); // geno
		data.addFile(eigenArgs.getPhenotypeFile(), eigenArgs.getSelectedPhenotype()); // pheno
		if (eigenArgs.isKeepFile())
			data.addFile(eigenArgs.getKeepFile()); // keep
		data.LineUpFiles();

		sf = new SampleFilter(pp.getPedigreeData(), data.getMatchSubjetList());
		pGM = new GenotypeMatrix(sf.getSample(), pp.getMapData(), cmdArgs);

		eigenGWASMT();
		printResult();
	}

	private void eigenGWASMT() {

		ChiSquaredDistributionImpl ci = new ChiSquaredDistributionImpl(1);

		int[] pIdx = this.data.getMatchedSubjectIdx(1);

		double[] Y = new double[pIdx.length];
		ArrayList<Double> pArray = NewIt.newArrayList();
		double threshold_ = 0.0D;

		for (int subjectIdx = 0; subjectIdx < Y.length; subjectIdx++) {
			Y[subjectIdx] = this.data.getVariable(1, pIdx[subjectIdx], this.traitIdx);
			threshold_ += Y[subjectIdx];
		}
		threshold_ /= Y.length;

		final double threshold = threshold_;

		int cpuNum = 1;
		if (eigenArgs.isThreadNum()) {
			cpuNum = eigenArgs.getThreadNum();
		}

		int markCnt = pGM.getNumMarker();
		double[] pVec = new double[markCnt];

		final int cpus = cpuNum < markCnt ? cpuNum : 1;

		Logger.printUserLog("Calculating eGWAS with " + cpus + (cpus == 1 ? " thread." : " threads."));

		final EigenGWASResult[] eRes = new EigenGWASResult[markCnt];
		Thread[] computeThreads = new Thread[cpus];
		final int[] taskProgresses = new int[cpus];
		final int taskSize = markCnt / cpus;

		for (int i = 0; i < cpus; ++i) {
			final int threadIndex = i;
			Thread thread = new Thread() {
				public void run() {

					int markStart = threadIndex * taskSize;
					int markEnd = threadIndex < (cpus - 1) ? (taskSize * (threadIndex+1)) : markCnt;

					int taskProgress = 0;

					for (int j = markStart; j < markEnd; j++) {

						SNP snp = pGM.getSNPList().get(j);

						SimpleRegression sReg = new SimpleRegression();
						double n1 = 0.0D;
						double n2 = 0.0D;
						double N = 0.0D;
						double freq1 = 0.0D;
						double freq2 = 0.0D;
						double freq = 0.0D;

						for (int k = 0; k < pGM.getNumIndivdial(); k++) {

							// int g = pGM.getAdditiveScoreOnFirstAllele(gIdx[j], i);
							int g = pGM.getAdditiveScoreOnFirstAllele(k, j);
							if (g != ConstValues.MISSING_GENOTYPE) {
								sReg.addData(g, Y[k]);
								if (Y[k] < threshold) {
									n1 += 1.0D;
									freq1 += g;
								} else {
									n2 += 1.0D;
									freq2 += g;
								}
								N += 1.0D;
								freq += g;
							}
						}

						freq1 /= 2.0D * n1;
						freq2 /= 2.0D * n2;
						freq /= 2.0D * N;

						double b, b_se;
						boolean isGood;
						if (freq == 0 || freq == 1 || (N < pGM.getNumIndivdial() * eigenArgs.getGENO())
								|| pGM.getAlleleVar(j) == 0) {
							b = Double.NaN;
							b_se = Double.NaN;
							isGood = false;
						} else {
							b = sReg.getSlope();
							b_se = sReg.getSlopeStdErr();
							isGood = true;
						}
						EigenGWASResult e1 = new EigenGWASResult(snp, freq, b, b_se, n1, freq1, n2, freq2, isGood);
						eRes[j] = e1;
//						pArray.add(e1.GetP());
						pVec[j] = e1.GetP();
						taskProgresses[threadIndex] = ++taskProgress;
					}
				}
			};
			thread.start();
			computeThreads[i] = thread;
		}

		Thread progressDisplayThread = new Thread() {
			public void run() {
				int totalProgress;
				do {
					totalProgress = IntStream.of(taskProgresses).sum();
					float percentage = Math.min(100f, (float) totalProgress / pGM.getNumMarker() * 100f);
					System.out.print(String.format(
							"\r[INFO] Calculating eGWAS, %.2f%% completed...",
							percentage));
					try {
						Thread.sleep(1000);
					} catch (InterruptedException e) {
					}
				} while (totalProgress < pGM.getNumMarker());
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

		for(int i = 0; i < pVec.length; i++) pArray.add(pVec[i]);
		eGWASResult = eRes;

		Collections.sort(pArray);
		int idx = (int) Math.ceil(pArray.size() / 2);

		if (Math.log10(pArray.size()) <= 6) {
			GC = new double[(int) (Math.log10(pArray.size()))];
			pQuantile = new double[(int) (Math.log10(pArray.size()))];
		} else {
			GC = new double[7];
			pQuantile = new double[7];
		}
		Logger.printUserLog("Median of p values is " + pArray.get(idx));

		try {
			double chisq = ci.inverseCumulativeProbability(1 - pArray.get(idx).doubleValue());
			lambdaGC = chisq / 0.4549;
			GC[0] = lambdaGC;
			pQuantile[0] = 0.5;
			for (int i = 1; i < GC.length; i++) {
				chisq = ci.inverseCumulativeProbability(Math.pow(10, -1 * i));
				GC[i] = ci.inverseCumulativeProbability(1 - pArray.get((int) (pArray.size() * Math.pow(10, -1 * i))));
				pQuantile[i] = Math.pow(10, -1 * i);
			}
		} catch (MathException e) {
			e.printStackTrace();
		}

		Logger.printUserLog("Lambda GC is: " + lambdaGC);

		if (eigenArgs.isGUI()) {
			PrintStream gui_file = null;
			gui_file = FileUtil.CreatePrintStream(eigenArgs.getOutRoot() + ".egwas.gui");
			gui_file.println(lambdaGC);
			gui_file.close();
		}
	}

	private void eigenGWAS() {

		ChiSquaredDistributionImpl ci = new ChiSquaredDistributionImpl(1);

		// int[] gIdx = this.data.getMatchedSubjectIdx(0);
		int[] pIdx = this.data.getMatchedSubjectIdx(1);

		double[] Y = new double[pIdx.length];
		ArrayList<Double> pArray = NewIt.newArrayList();
		double threshold = 0.0D;

		for (int subjectIdx = 0; subjectIdx < Y.length; subjectIdx++) {
			Y[subjectIdx] = this.data.getVariable(1, pIdx[subjectIdx], this.traitIdx);
			threshold += Y[subjectIdx];
		}
		threshold /= Y.length;

		for (int i = 0; i < pGM.getNumMarker(); i++) {
			SNP snp = pGM.getSNPList().get(i);

			SimpleRegression sReg = new SimpleRegression();
			double n1 = 0.0D;
			double n2 = 0.0D;
			double N = 0.0D;
			double freq1 = 0.0D;
			double freq2 = 0.0D;
			double freq = 0.0D;

			for (int j = 0; j < pGM.getNumIndivdial(); j++) {

				// int g = pGM.getAdditiveScoreOnFirstAllele(gIdx[j], i);
				int g = pGM.getAdditiveScoreOnFirstAllele(j, i);
				if (g != ConstValues.MISSING_GENOTYPE) {
					sReg.addData(g, Y[j]);
					if (Y[j] < threshold) {
						n1 += 1.0D;
						freq1 += g;
					} else {
						n2 += 1.0D;
						freq2 += g;
					}
					N += 1.0D;
					freq += g;
				}
			}

			freq1 /= 2.0D * n1;
			freq2 /= 2.0D * n2;
			freq /= 2.0D * N;

			double b, b_se;
			boolean isGood;
			if (freq == 0 || freq == 1 || (N < pGM.getNumIndivdial() * eigenArgs.getGENO())
					|| pGM.getAlleleVar(i) == 0) {
				b = Double.NaN;
				b_se = Double.NaN;
				isGood = false;
			} else {
				b = sReg.getSlope();
				b_se = sReg.getSlopeStdErr();
				isGood = true;
			}
			EigenGWASResult e1 = new EigenGWASResult(snp, freq, b, b_se, n1, freq1, n2, freq2, isGood);
			eGWASResult[i] = e1;
			pArray.add(e1.GetP());
		}
		Collections.sort(pArray);
		int idx = (int) Math.ceil(pArray.size() / 2);

		if (Math.log10(pArray.size()) <= 6) {
			GC = new double[(int) (Math.log10(pArray.size()))];
			pQuantile = new double[(int) (Math.log10(pArray.size()))];
		} else {
			GC = new double[7];
			pQuantile = new double[7];
		}
		Logger.printUserLog("Median of p values is " + pArray.get(idx));

		try {
			double chisq = ci.inverseCumulativeProbability(1 - pArray.get(idx).doubleValue());
			lambdaGC = chisq / 0.4549;
			GC[0] = lambdaGC;
			pQuantile[0] = 0.5;
			for (int i = 1; i < GC.length; i++) {
				chisq = ci.inverseCumulativeProbability(Math.pow(10, -1 * i));
				GC[i] = ci.inverseCumulativeProbability(1 - pArray.get((int) (pArray.size() * Math.pow(10, -1 * i))));
				pQuantile[i] = Math.pow(10, -1 * i);
			}
		} catch (MathException e) {
			e.printStackTrace();
		}

		Logger.printUserLog("Lambda GC is: " + lambdaGC);

		if (eigenArgs.isGUI()) {
			PrintStream gui_file = null;
			gui_file = FileUtil.CreatePrintStream(eigenArgs.getOutRoot() + ".egwas.gui");
			gui_file.println(lambdaGC);
			gui_file.close();
		}
	}

	public void printResult() {

		PrintStream eGWAS = FileUtil.CreatePrintStream(this.eigenArgs.getOutRoot() + ".egwas");
		if (eigenArgs.isTAB()) {
			eGWAS.println("SNP\tRefAllele\tBeta");
		} else {
			eGWAS.println("SNP\tCHR\tBP\tRefAllele\tAltAllele\tfreq\tBeta\tSE\tChi\tP\tPGC\tn1\tfreq1\tn2\tfreq2\tFst");
		}

		for (int i = 0; i < eGWASResult.length; i++) {
			EigenGWASResult e1 = eGWASResult[i];
			int j = 0;
			for (; j < GC.length-1; j++) {
				if (e1.GetP() >= pQuantile[j]) {
					break;
				}
			}
			if (eigenArgs.isTAB()) {
				eGWAS.println(e1.printEGWASTab(GC[j]));
			} else {
				eGWAS.println(e1.printEGWASResult(GC[j]));
				//eGWAS.println(e1.printEGWASResult(lambdaGC));
			}
		}
		eGWAS.close();
	}
}