package gear.subcommands.eigengwasdom;

import gear.data.InputDataSet2;
import gear.family.pedigree.file.SNP;
import gear.family.plink.PLINKParser;
import gear.family.popstat.GenotypeMatrix;
import gear.family.qc.rowqc.SampleFilter;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;
import gear.util.pop.PopStat;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ChiSquaredDistributionImpl;
import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;

public class EigenGWASDomImpl extends CommandImpl {
	private EigenGWASDomCommandArguments eigenArgs;
	private SampleFilter sf;
	private GenotypeMatrix pGM;
	private int traitIdx;
	private InputDataSet2 data = new InputDataSet2();
	private ArrayList<EigenGWASDomResult> eGWASResult = NewIt.newArrayList();

	private double lambdaGC = 1;
	private double lambdaGCDom = 1;
	private int monoLoci = 0;
	private int singularLoci = 0;

	public void execute(CommandArguments cmdArgs) {
		this.eigenArgs = ((EigenGWASDomCommandArguments) cmdArgs);

		this.traitIdx = this.eigenArgs.getMpheno()[0];
		data.addFile(this.eigenArgs.getFam());
		data.addFile(this.eigenArgs.getPhenotypeFile(), this.eigenArgs.getMpheno());
		if (eigenArgs.isKeepFile()) data.addFile(this.eigenArgs.getKeepFile());
		data.LineUpFiles();

		PLINKParser pp = PLINKParser.parse(this.eigenArgs);
		sf = new SampleFilter(pp.getPedigreeData(), cmdArgs);
		pGM = new GenotypeMatrix(sf.getSample(), pp.getMapData(), cmdArgs);

		eigenGWAS();
		printResult();
	}

	private void eigenGWAS() {
		ChiSquaredDistributionImpl ci = new ChiSquaredDistributionImpl(1);

		int[] gIdx = this.data.getMatchedSubjectIdx(0);
		int[] pIdx = this.data.getMatchedSubjectIdx(1);

		double[] Y = new double[pIdx.length];

		ArrayList<Double> pArray = NewIt.newArrayList();
		ArrayList<Double> pDomArray = NewIt.newArrayList();

		double threshold = 0.0D;

		for (int subjectIdx = 0; subjectIdx < Y.length; subjectIdx++) {
			Y[subjectIdx] = this.data.getVariable(1, pIdx[subjectIdx], this.traitIdx);
			threshold += Y[subjectIdx];
		}
		threshold /= Y.length;

		// PrintStream eGWAS =
		// FileUtil.CreatePrintStream(this.eigenArgs.getOutRoot() + ".egwas");
		// eGWAS.println("SNP\tCHR\tBP\tRefAllele\tAltAllele\tfreq\tBeta\tSE\tChi\tP\tPGC\tn1\tfreq1\tn2\tfreq2\tFst");

		double[][] gfreq = PopStat.calGenoFrequency(pGM);
		for (int i = 0; i < pGM.getNumMarker(); i++) {
			SNP snp = pGM.getSNPList().get(i);

			double[][] x = new double[gIdx.length][2];
			ArrayList<Integer> IDX = NewIt.newArrayList();
			double n1 = 0.0D;
			double n2 = 0.0D;
			double N = 0.0D;
			double freq1 = 0.0D;
			double freq2 = 0.0D;
			double freq = 0.0D;
			for (int j = 0; j < gIdx.length; j++) {
				int g = pGM.getAdditiveScoreOnFirstAllele(gIdx[j], i);
				if (g != 3) {
					IDX.add(j);
					if (g == 0) {
						/*
						 * x[j][0] = -1 * (2 * gfreq[i][2] + gfreq[i][1]); x[j][1] = 1/ (8 *
						 * gfreq[i][0]);
						 */
						x[j][0] = -gfreq[i][1] - 2 * gfreq[i][2];
						x[j][1] = (2 * gfreq[i][1] * gfreq[i][2]) / (gfreq[i][0] + gfreq[i][2]
								- (gfreq[i][0] - gfreq[i][2]) * (gfreq[i][0] - gfreq[i][2]));
					} else if (g == 1) {
						/*
						 * x[j][0] = gfreq[i][0] - gfreq[i][2]; x[j][1] = -1 / (4 * gfreq[i][1]);
						 */
						x[j][0] = 1 - gfreq[i][1] - 2 * gfreq[i][2];
						x[j][1] = -1 * (4 * gfreq[i][0] * gfreq[i][2]) / (gfreq[i][0] + gfreq[i][2]
								- (gfreq[i][0] - gfreq[i][2]) * (gfreq[i][0] - gfreq[i][2]));
					} else {
						/*
						 * x[j][0] = -1 * (2 - 1*gfreq[i][1] - 2*gfreq[i][2]); x[j][1] = -1 / (8 *
						 * gfreq[i][2]);
						 */
						x[j][0] = 2 - gfreq[i][1] - 2 * gfreq[i][2];
						x[j][1] = (2 * gfreq[i][0] * gfreq[i][1]) / (gfreq[i][0] + gfreq[i][2]
								- (gfreq[i][0] - gfreq[i][2]) * (gfreq[i][0] - gfreq[i][2]));
					}

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

				freq1 /= 2.0D * n1;
				freq2 /= 2.0D * n2;
				freq /= 2.0D * N;

				double fst = 2 * (n1 / N * (freq1 - freq) * (freq1 - freq) + n2 / N * (freq2 - freq) * (freq2 - freq))
						/ (freq * (1.0D - freq));

				double[][] x1 = new double[IDX.size()][3];
				double[] y1 = new double[IDX.size()];
				for (int k1 = 0; k1 < x1.length; k1++) {
					x1[k1][0] = 1;
					x1[k1][1] = x[IDX.get(k1)][0];
					x1[k1][2] = x[IDX.get(k1)][1];
					y1[k1] = Y[IDX.get(k1)];
				}

				RealMatrix g1 = new Array2DRowRealMatrix(x1);
				RealMatrix XtX = g1.transpose().multiply(g1);

				boolean isNonSingular = (new LUDecompositionImpl(XtX)).getSolver().isNonSingular();

				if (!isNonSingular) {
					Logger.printUserLog("Model is singular for '" + snp.getName() + "'; skipped.");
					singularLoci++;
					continue;
				} else {
					RealMatrix XtX_inv = (new LUDecompositionImpl(XtX)).getSolver().getInverse();

					RealMatrix g1_tran = g1.transpose();
					RealMatrix y = new Array2DRowRealMatrix(y1);
					RealMatrix B = XtX_inv.multiply(g1_tran).multiply(y);
					RealMatrix SST = y.transpose().multiply(y);
					RealMatrix SSR = B.transpose().multiply(g1_tran).multiply(y);
					double sse = (SST.getEntry(0, 0) - SSR.getEntry(0, 0))
							/ (y.getRowDimension() - B.getRowDimension());

					RealMatrix BV = XtX_inv.scalarMultiply(sse);
					if (BV.getEntry(1, 1) > 0 && BV.getEntry(2, 2) > 0) {
						EigenGWASDomResult e1 = new EigenGWASDomResult(snp, freq, B.getEntry(1, 0),
								Math.sqrt(BV.getEntry(1, 1)), B.getEntry(2, 0), Math.sqrt(BV.getEntry(2, 2)), n1, freq1,
								n2, freq2, fst);
						eGWASResult.add(e1);
						pArray.add(e1.GetP());
						pDomArray.add(e1.GetDomP());
					}
				}
			}
		}
		Collections.sort(pArray);
		int idx = (int) Math.ceil(pArray.size() / 2);

		Collections.sort(pDomArray);
		int idxDom = (int) Math.ceil(pDomArray.size() / 2);

		if (monoLoci > 1) {
			Logger.printUserLog("Removed " + monoLoci + " monomorphic loci.");
		} else if (monoLoci == 1) {
			Logger.printUserLog("Removed " + monoLoci + " monomorphic locus.");
		}

		if (singularLoci > 1) {
			Logger.printUserLog("Removed " + singularLoci + " singular loci.");
		} else if (singularLoci == 1) {
			Logger.printUserLog("Removed " + singularLoci + " singular locus.");
		}
		Logger.printUserLog("Median of p values (additive) is " + pArray.get(idx));
		Logger.printUserLog("Median of p values (dominance) is " + pDomArray.get(idx));

		try {
			double chisq = ci.inverseCumulativeProbability(1 - pArray.get(idx).doubleValue());
			lambdaGC = chisq / 0.4549;

			double chisqDom = ci.inverseCumulativeProbability(1 - pDomArray.get(idxDom).doubleValue());
			lambdaGCDom = chisqDom / 0.4549;
		} catch (MathException e) {
			e.printStackTrace();
		}

		Logger.printUserLog("Lambda GC (additive) is: " + lambdaGC);
		Logger.printUserLog("Lambda GC (dominance) is: " + lambdaGCDom);
	}

	public void printResult() {
		PrintStream eGWAS = FileUtil.CreatePrintStream(this.eigenArgs.getOutRoot() + ".egwasd");
		eGWAS.println(
				"SNP\tCHR\tBP\tRefAllele\tAltAllele\tfreq\tBeta\tSE\tChi\tP\tPGC\tDom\tDomSE\tpDom\tn1\tfreq1\tn2\tfreq2\tFst");

		for (int i = 0; i < eGWASResult.size(); i++) {
			EigenGWASDomResult e1 = eGWASResult.get(i);
			eGWAS.println(e1.printEGWASResult(lambdaGC));
		}
		eGWAS.close();
	}
}