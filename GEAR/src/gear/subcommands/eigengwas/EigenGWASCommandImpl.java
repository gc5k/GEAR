package gear.subcommands.eigengwas;

import gear.data.InputDataSet2;
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
import gear.util.NewIt;

import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ChiSquaredDistributionImpl;
import org.apache.commons.math.stat.regression.SimpleRegression;

public class EigenGWASCommandImpl extends CommandImpl {
	private EigenGWASCommandArguments eigenArgs;
	private SampleFilter sf;

	private MapFile map;
	private SNPFilter snpFilter;
	private BEDReader bed;

	private int traitIdx;
	private InputDataSet2 data = new InputDataSet2();
	private ArrayList<EigenGWASResult> eGWASResult = NewIt.newArrayList();

	private double lambdaGC = 1;
	private double[] GC;
	private double[] pQuantile;
	private DecimalFormat fmt1 = new DecimalFormat("0.0000");
	private DecimalFormat fmt2 = new DecimalFormat("0.00E000");

	public void execute(CommandArguments cmdArgs) {

		eigenArgs = (EigenGWASCommandArguments) cmdArgs;

		PLINKParser pp = PLINKParser.create(eigenArgs);
		pp.parseSmallFiles();
		map = pp.getMapData();
		snpFilter = pp.getSNPFilter();

		traitIdx = eigenArgs.getSelectedPhenotype(0);
		data.addFile(eigenArgs.getFam()); // geno
		data.addFile(eigenArgs.getPhenotypeFile(), eigenArgs.getSelectedPhenotype()); // pheno
		if (eigenArgs.isKeepFile())
			data.addFile(eigenArgs.getKeepFile()); // keep
		data.LineUpFiles();

		sf = new SampleFilter(pp.getPedigreeData(), data.getMatchSubjetList());
		bed = (BEDReader) pp.getPedigreeData();
		eigenGWASBed();
		printResult();
	}

	private void eigenGWASBed() {
		ChiSquaredDistributionImpl ci = new ChiSquaredDistributionImpl(1);
		SNPFilterPostQC snpPostQC = new SNPFilterPostQC(eigenArgs);
		ArrayList<Hukou> hkBook = bed.getHukouBook();

		double threshold = 0.0D;
		int[] pIdx = this.data.getMatchedSubjectIdx(1);
		double[] Y = new double[pIdx.length];
		ArrayList<Double> pArray = NewIt.newArrayList();

		for (int subjectIdx = 0; subjectIdx < Y.length; subjectIdx++) {
			Y[subjectIdx] = this.data.getVariable(1, pIdx[subjectIdx], this.traitIdx);
			threshold += Y[subjectIdx];
		}
		threshold /= Y.length;

		int numMarkers = map.getMarkerNumberOriginal();
		int numSamples = bed.getNumIndividuals();
		int workingSnpIndex = 0;


		for (int i = 0; i < numMarkers; ++i) {

			if (snpFilter.isSnpIncluded(i)) {

				SimpleRegression sReg = new SimpleRegression();
				int missingCnt = 0;
				int sum = 0;
				int squareSum = 0;

				int[] FmissingCnt = { 0, 0 };
				int[] Fsum = { 0, 0 };
				int[] Fn = { 0, 0 };

				for (int j = 0; j < numSamples; j += 4) {
					int nextByte = bed.readNextByte();
					int indCnt = j;
					for (int k = 0; k < 8; k += 2) {
						if (indCnt == numSamples)
							break;
						if (!hkBook.get(indCnt++).isAvailable()) {
							continue;
						}

						boolean isMissing = false;
						int gIdx = Y[indCnt - 1] < threshold ? 0 : 1;

						int genotype = (nextByte >> k) & 0b11;

						switch (genotype) {
						case PLINKBinaryParser.HOMOZYGOTE_FIRST:
							break;
						case PLINKBinaryParser.HETEROZYGOTE:
							sum += 1;
							squareSum += 1;
							Fsum[gIdx] += 1;
							break;
						case PLINKBinaryParser.HOMOZYGOTE_SECOND:
							sum += 2;
							squareSum += 4;
							Fsum[gIdx] += 2;
							break;
						case PLINKBinaryParser.MISSING_GENOTYPE:
							++missingCnt;
							isMissing = true;
							++FmissingCnt[gIdx];
							break;
						}
						if (!isMissing) {
							Fn[gIdx]++;
							sReg.addData(genotype, Y[indCnt - 1]);
						}
					}
				}
				int validSampleCnt = numSamples - missingCnt;
				double variance = 0;
				double average = 0;
				if (validSampleCnt > 2) {
					average = (double) sum / (2.0*validSampleCnt);
					variance = (squareSum - validSampleCnt * average * average) / (validSampleCnt - 1);
				}

				double freq = 1 - average;
				double maf = freq < 0.5 ? freq : (1 - freq);
				double freq1 = 1 - (Fsum[0] * 1.0D) / (2.0D * Fn[0]);
				double freq2 = 1 - (Fsum[1] * 1.0D) / (2.0D * Fn[1]);

				boolean isPassPostQC = snpPostQC.isPassPostQC(maf, missingCnt * 1.0D / numSamples);
				if (!isPassPostQC) continue;

				double b, b_se;
				boolean isGood = true;
				if (freq == 0 || variance == 0) {
					b = Double.NaN;
					b_se = Double.NaN;
					isGood = false;
				} else {
					b = sReg.getSlope();
					b_se = sReg.getSlopeStdErr();
				}
				SNP snp = map.getSNP(workingSnpIndex++);
				EigenGWASResult e1 = new EigenGWASResult(snp, freq, b, b_se, Fn[0], freq1, Fn[1], freq2, isGood);
				eGWASResult.add(e1);
				pArray.add(e1.GetP());
			}
		}
		snpPostQC.printPostQCSummary();
		Collections.sort(pArray);
		int idx = (int) Math.ceil(pArray.size() / 2);

		if (Math.log10(pArray.size()) <= 6) {
			GC = new double[(int) (Math.log10(pArray.size()))];
			pQuantile = new double[(int) (Math.log10(pArray.size()))];
		} else {
			GC = new double[7];
			pQuantile = new double[7];
		}
		double pMedian = pArray.get(idx);
		Logger.printUserLog("Median of p values is " + (pMedian > 0.0001? fmt1.format(pMedian):fmt2.format(pMedian))+ " for " + pArray.size() +" analyzed SNPs.");
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

		Logger.printUserLog("Lambda GC is: " + fmt1.format(lambdaGC));

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

		for (int i = 0; i < eGWASResult.size(); i++) {
			EigenGWASResult e1 = eGWASResult.get(i);
			int j = 0;
			for (; j < GC.length; j++) {
				if (e1.GetP() >= pQuantile[j]) {
					break;
				}
			}
			if (eigenArgs.isTAB()) {
				// eGWAS.println(e1.printEGWASTab(lambdaGC));
				eGWAS.println(e1.printEGWASTab(GC[j]));
			} else {
				// eGWAS.println(e1.printEGWASResult(GC[j]));
				eGWAS.println(e1.printEGWASResult(lambdaGC));
			}
		}
		eGWAS.close();
	}
}