package gear.subcommands.eigengwas;

import gear.data.InputDataSet2;
import gear.family.GenoMatrix.GenotypeMatrix;
import gear.family.pedigree.file.SNP;
import gear.family.plink.PLINKParser;
import gear.family.qc.rowqc.SampleFilter;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ChiSquaredDistributionImpl;
import org.apache.commons.math.stat.regression.SimpleRegression;

public class EigenGWASCommandImpl extends CommandImpl {
	private EigenGWASCommandArguments eigenArgs;
	private SampleFilter sf;
	private GenotypeMatrix pGM;

	private int traitIdx;
	private InputDataSet2 data = new InputDataSet2();
	private ArrayList<EigenGWASResult> eGWASResult = NewIt.newArrayList();

	private double lambdaGC = 1;

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

		eigenGWAS();
		printResult();
	}

	private void eigenGWAS() {
		// ////TEST
		// File f = new File(this.eigenArgs.getOutRoot() + ".egwas");
		// BufferedWriter bw = null;
		// OutputStreamWriter write = null;
		// try {
		// f.createNewFile();
		// write = new OutputStreamWriter(new FileOutputStream(f));
		// bw = new BufferedWriter(write);
		// bw.write("SNP\tCHR\tBP\tRefAllele\tAltAllele\tfreq\tBeta\tSE\tChi\tP\tPGC\tn1\tfreq1\tn2\tfreq2\tFst\n");
		// } catch (IOException e2) {
		// e2.printStackTrace();
		// }
		// ////TEST

		ChiSquaredDistributionImpl ci = new ChiSquaredDistributionImpl(1);

//		int[] gIdx = this.data.getMatchedSubjectIdx(0);
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
				
//				int g = pGM.getAdditiveScoreOnFirstAllele(gIdx[j], i);
				int g = pGM.getAdditiveScoreOnFirstAllele(j, i);
				if (g != 3) {
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

			double fst, b, b_se;
			boolean isGood;
			if (freq == 0 || freq == 1 || (N < pGM.getNumIndivdial() * eigenArgs.getGENO())
					|| pGM.getAlleleVar(i) == 0) {
				fst = Double.NaN;
				b = Double.NaN;
				b_se = Double.NaN;
				isGood = false;
				// continue;
			} else {
				fst = 2 * (n1 / N * (freq1 - freq) * (freq1 - freq) + n2 / N * (freq2 - freq) * (freq2 - freq))
						/ (freq * (1.0D - freq));
				b = sReg.getSlope();
				b_se = sReg.getSlopeStdErr();
				isGood = true;
			}
			EigenGWASResult e1 = new EigenGWASResult(snp, freq, b, b_se, n1, freq1, n2, freq2, fst, isGood);
			eGWASResult.add(e1);
			pArray.add(e1.GetP());

			// ///TEST
			// try {
			// bw.write(e1.printEGWASResult(1) + "\n");
			// if (i > 10000)
			// {
			// bw.flush();
			// }
			//
			// } catch (IOException e) {
			// e.printStackTrace();
			// }
			// ///TEST
		}
		Collections.sort(pArray);
		int idx = (int) Math.ceil(pArray.size() / 2);

		Logger.printUserLog("Median of p values is " + pArray.get(idx));

		try {
			double chisq = ci.inverseCumulativeProbability(1 - pArray.get(idx).doubleValue());
			lambdaGC = chisq / 0.4549;
		} catch (MathException e) {
			e.printStackTrace();
		}

		Logger.printUserLog("Lambda GC is: " + lambdaGC);

		// ////TEST
		// try {
		// write.close();
		// bw.close();
		//
		// } catch (IOException e) {
		// e.printStackTrace();
		// }
		// ////TEST
	}

	public void printResult() {

		PrintStream eGWAS = FileUtil.CreatePrintStream(this.eigenArgs.getOutRoot() + ".egwas");
		eGWAS.println("SNP\tCHR\tBP\tRefAllele\tAltAllele\tfreq\tBeta\tSE\tChi\tP\tPGC\tn1\tfreq1\tn2\tfreq2\tFst");

		for (int i = 0; i < eGWASResult.size(); i++) {
			EigenGWASResult e1 = eGWASResult.get(i);
			eGWAS.println(e1.printEGWASResult(lambdaGC));
		}
		eGWAS.close();
	}
}