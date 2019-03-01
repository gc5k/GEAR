package gear.subcommands.eigengwas;

import gear.ConstValues;
import gear.data.InputDataSet2;
import gear.family.GenoMatrix.GenotypeMatrix;
import gear.family.pedigree.Hukou;
import gear.family.pedigree.file.BEDReader;
import gear.family.pedigree.file.MapFile;
import gear.family.pedigree.file.SNP;
import gear.family.plink.PLINKBinaryParser;
import gear.family.plink.PLINKParser;
import gear.family.qc.colqc.SNPFilter;
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

	private MapFile map;
	private SNPFilter snpFilter;
	private BEDReader bed;
	
	private int traitIdx;
	private InputDataSet2 data = new InputDataSet2();
	private ArrayList<EigenGWASResult> eGWASResult = NewIt.newArrayList();

	private double lambdaGC = 1;
	private double[] GC;
	private double[] pQuantile;

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
		pGM = new GenotypeMatrix(sf.getSample(), pp.getMapData(), cmdArgs);

//		eigenGWAS();
		if (pp.getPedigreeData() instanceof BEDReader) {
			bed = (BEDReader) pp.getPedigreeData();
			eigenGWASBed();
			printResult();
		}
	}

	private void eigenGWASBed() {
		ChiSquaredDistributionImpl ci = new ChiSquaredDistributionImpl(1);

		int[] pIdx = this.data.getMatchedSubjectIdx(1);
		double[] Y = new double[pIdx.length];
		ArrayList<Double> pArray = NewIt.newArrayList();
		double threshold = 0.0D;

		for (int subjectIdx = 0; subjectIdx < Y.length; subjectIdx++) {
			Y[subjectIdx] = this.data.getVariable(1, pIdx[subjectIdx], this.traitIdx);
			threshold += Y[subjectIdx];
		}
		threshold /= Y.length;
		
		int numMarkers = map.getMarkerNumberOriginal();
		int numSamples = bed.getNumIndividuals();
		int workingSnpIndex = 0;

		ArrayList<Hukou> hkBook = bed.getHukouBook();

		for (int i = 0; i < numMarkers; ++i) {

			if (snpFilter.isSnpIncluded(i)) {

				SimpleRegression sReg = new SimpleRegression();
				int missingCnt = 0;
				int sum = 0;
				int squareSum = 0;

				int[] FmissingCnt = {0, 0};
				int[] Fsum = {0, 0};
				int[] Fn = {0, 0};
				
				for (int j = 0; j < numSamples; j += 4) {
					int nextByte = bed.readNextByte();
					int indCnt = j;
					for (int k = 0; k < 8; k += 2) {
						if(indCnt == numSamples) break;
						boolean isMissing = false;
						if (!hkBook.get(indCnt++).isAvailable()) {
							continue;
						}
						int gIdx = Y[indCnt-1] < threshold ? 0:1;

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
							sReg.addData(genotype, Y[indCnt-1]); 
						}
					}
				}
				int validSampleCnt = numSamples - missingCnt;
				double variance = 0;
				if (validSampleCnt > 2) {
					double average = (double)sum / validSampleCnt;
					variance = (squareSum - validSampleCnt * average * average) / (validSampleCnt - 1);
				}

				double freq = 1-(sum*1.0D)/(2*validSampleCnt);
				double freq1 = 1-(Fsum[0]*1.0D)/(2.0D*Fn[0]);
				double freq2 = 1-(Fsum[1]*1.0D)/(2.0D*Fn[1]);

				double b, b_se;
				boolean isGood;
				if (freq == 0 || freq == 1 || ((Fn[0]+Fn[1]) < numSamples * eigenArgs.getGENO())
					|| variance == 0) {
					b = Double.NaN;
					b_se = Double.NaN;
					isGood = false;
				} else {
					b = sReg.getSlope();
					b_se = sReg.getSlopeStdErr();
					isGood = true;
				}
//				System.out.println(snp.getName() + " freq= "+ freq + " b="+b + " b_se=" + b_se + " f1="+ freq1 + " f2=" + freq2);
				SNP snp = map.getSNP(workingSnpIndex++);

				EigenGWASResult e1 = new EigenGWASResult(snp, freq, b, b_se, Fn[0], freq1, Fn[1], freq2, isGood);
				eGWASResult.add(e1);
				pArray.add(e1.GetP());
			}
		}
		Collections.sort(pArray);
		int idx = (int) Math.ceil(pArray.size() / 2);

		if(Math.log10(pArray.size()) <= 6) {
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
			for(int i = 1; i < GC.length; i++) {
				chisq = ci.inverseCumulativeProbability(Math.pow(10, -1*i));
				GC[i] = ci.inverseCumulativeProbability(1 - pArray.get((int) (pArray.size()*Math.pow(10, -1*i))) );
				pQuantile[i] = Math.pow(10, -1*i);
			}
		} catch (MathException e) {
			e.printStackTrace();
		}

		Logger.printUserLog("Lambda GC is: " + lambdaGC);

		if (eigenArgs.isGUI()) {
			PrintStream gui_file = null;
			gui_file = FileUtil.CreatePrintStream(eigenArgs.getOutRoot()+".egwas.gui");
			gui_file.println(lambdaGC);
			gui_file.close();
		}
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

//	private void eigenGWAS() {
//		// ////TEST
//		// File f = new File(this.eigenArgs.getOutRoot() + ".egwas");
//		// BufferedWriter bw = null;
//		// OutputStreamWriter write = null;
//		// try {
//		// f.createNewFile();
//		// write = new OutputStreamWriter(new FileOutputStream(f));
//		// bw = new BufferedWriter(write);
//		// bw.write("SNP\tCHR\tBP\tRefAllele\tAltAllele\tfreq\tBeta\tSE\tChi\tP\tPGC\tn1\tfreq1\tn2\tfreq2\tFst\n");
//		// } catch (IOException e2) {
//		// e2.printStackTrace();
//		// }
//		// ////TEST
//
//		ChiSquaredDistributionImpl ci = new ChiSquaredDistributionImpl(1);
//
////		int[] gIdx = this.data.getMatchedSubjectIdx(0);
//		int[] pIdx = this.data.getMatchedSubjectIdx(1);
//
//		double[] Y = new double[pIdx.length];
//		ArrayList<Double> pArray = NewIt.newArrayList();
//		double threshold = 0.0;
//
//		for (int subjectIdx = 0; subjectIdx < Y.length; subjectIdx++) {
//			Y[subjectIdx] = this.data.getVariable(1, pIdx[subjectIdx], this.traitIdx);
//			threshold += Y[subjectIdx];
//		}
//		threshold /= Y.length;
//
//		for (int i = 0; i < pGM.getNumMarker(); i++) {
//			SNP snp = pGM.getSNPList().get(i);
//
//			SimpleRegression sReg = new SimpleRegression();
//			double n1 = 0.0D;
//			double n2 = 0.0D;
//			double N = 0.0D;
//			double freq1 = 0.0D;
//			double freq2 = 0.0D;
//			double freq = 0.0D;
//
//			for (int j = 0; j < pGM.getNumIndivdial(); j++) {
//				
////				int g = pGM.getAdditiveScoreOnFirstAllele(gIdx[j], i);
//				int g = pGM.getAdditiveScoreOnFirstAllele(j, i);
//				if (g != ConstValues.MISSING_GENOTYPE) {
//					sReg.addData(g, Y[j]);
//					if (Y[j] < threshold) {
//						n1 += 1.0D;
//						freq1 += g;
//					} else {
//						n2 += 1.0D;
//						freq2 += g;
//					}
//					N += 1.0D;
//					freq += g;
//				}
//			}
//
//			freq1 /= 2.0D * n1;
//			freq2 /= 2.0D * n2;
//			freq /= 2.0D * N;
//
//			double b, b_se;
//			boolean isGood;
//			if (freq == 0 || freq == 1 || (N < pGM.getNumIndivdial() * eigenArgs.getGENO())
//					|| pGM.getAlleleVar(i) == 0) {
//				b = Double.NaN;
//				b_se = Double.NaN;
//				isGood = false;
//			} else {
//				b = sReg.getSlope();
//				b_se = sReg.getSlopeStdErr();
//				isGood = true;
//			}
//			EigenGWASResult e1 = new EigenGWASResult(snp, freq, b, b_se, n1, freq1, n2, freq2, isGood);
//			eGWASResult.add(e1);
//			pArray.add(e1.GetP());
//
//			// ///TEST
//			// try {
//			// bw.write(e1.printEGWASResult(1) + "\n");
//			// if (i > 10000)
//			// {
//			// bw.flush();
//			// }
//			//
//			// } catch (IOException e) {
//			// e.printStackTrace();
//			// }
//			// ///TEST
//		}
//		Collections.sort(pArray);
//		int idx = (int) Math.ceil(pArray.size() / 2);
//
//		if(Math.log10(pArray.size()) <= 6) {
//			GC = new double[(int) (Math.log10(pArray.size()))];
//			pQuantile = new double[(int) (Math.log10(pArray.size()))];
//		} else {
//			GC = new double[7];
//			pQuantile = new double[7];
//		}
//		Logger.printUserLog("Median of p values is " + pArray.get(idx));
//
//		try {
//			double chisq = ci.inverseCumulativeProbability(1 - pArray.get(idx).doubleValue());
//			lambdaGC = chisq / 0.4549;
//			GC[0] = lambdaGC;
//			pQuantile[0] = 0.5;
//			for(int i = 1; i < GC.length; i++) {
//				chisq = ci.inverseCumulativeProbability(Math.pow(10, -1*i));
//				GC[i] = ci.inverseCumulativeProbability(1 - pArray.get((int) (pArray.size()*Math.pow(10, -1*i))) );
//				pQuantile[i] = Math.pow(10, -1*i);
//			}
//		} catch (MathException e) {
//			e.printStackTrace();
//		}
//
//		Logger.printUserLog("Lambda GC is: " + lambdaGC);
//
//		if (eigenArgs.isGUI()) {
//			PrintStream gui_file = null;
//			gui_file = FileUtil.CreatePrintStream(eigenArgs.getOutRoot()+".egwas.gui");
//			gui_file.println(lambdaGC);
//			gui_file.close();
//		}
//		// ////TEST
//		// try {
//		// write.close();
//		// bw.close();
//		//
//		// } catch (IOException e) {
//		// e.printStackTrace();
//		// }
//		// ////TEST
//	}

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
//				eGWAS.println(e1.printEGWASTab(lambdaGC));
				eGWAS.println(e1.printEGWASTab(GC[j]));
			} else {
//				eGWAS.println(e1.printEGWASResult(GC[j]));
				eGWAS.println(e1.printEGWASResult(lambdaGC));				
			}
		}
		eGWAS.close();
	}
}