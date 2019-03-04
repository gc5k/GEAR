package gear.subcommands.oath.nss;

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
import gear.subcommands.oath.OATHConst;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ChiSquaredDistributionImpl;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.stat.regression.SimpleRegression;
import org.apache.commons.math.stat.correlation.PearsonsCorrelation;

public class NSSCommandImpl extends CommandImpl {
	private NSSCommandArguments nssArgs;
	private SampleFilter sf;
	
	private MapFile map;
	private SNPFilter snpFilter;
	private BEDReader bed;

	private int[] traitIdx;
	private int[] covIdx;
	private int N;
	private InputDataSet2 data = null;
	private ArrayList< ArrayList<NSSGWASResult>> nssResultList; 
	private ArrayList<ArrayList<Double>> pArrayList;

	private double lambdaGC = 1;

//	private int famFileIdx = 0;
	private int pheFileIdx = 1;
	private int covFileIdx = 2;

	private DecimalFormat fmt1 = new DecimalFormat("0.0000");
	private DecimalFormat fmt2 = new DecimalFormat("0.00E000");

	public void execute(CommandArguments cmdArgs) {
		this.nssArgs = ((NSSCommandArguments) cmdArgs);

		this.traitIdx = this.nssArgs.getSelectedPhenotype();
		this.covIdx = this.nssArgs.getCovNumber();

		data = new InputDataSet2();
		data.addFile(this.nssArgs.getFam());		
		data.addFile(this.nssArgs.getPhenotypeFile(), this.nssArgs.getSelectedPhenotype());

		data.addFile(this.nssArgs.getCovFile(), this.nssArgs.getCovNumber());
		if (this.nssArgs.isKeepFile())
			data.addFile(this.nssArgs.getKeepFile());
		data.LineUpFiles();

		this.N = data.getNumberOfSubjects();

		PLINKParser pp = PLINKParser.create(nssArgs);
		pp.parseSmallFiles();
		map = pp.getMapData();
		snpFilter = pp.getSNPFilter();

		sf = new SampleFilter(pp.getPedigreeData(), data.getMatchSubjetList());
		bed = (BEDReader) pp.getPedigreeData();

		nssResultList = NewIt.newArrayList();
		pArrayList = NewIt.newArrayList();

		naiveGWASBed();
		printResult();

		String Fout = nssArgs.getOutRoot() + ".list.nss";
		PrintStream nssList = FileUtil.CreatePrintStream(Fout);
		for (int i = 0; i < traitIdx.length; i++) {
			nssList.println(nssArgs.getOutRoot() + ".p." + (traitIdx[i] + 1) + ".nss.gz");
		}
		for (int i = 0; i < covIdx.length; i++) {
			nssList.println(nssArgs.getOutRoot() + ".c." + (covIdx[i] + 1) + ".nss.gz");
		}

		nssList.close();
		printCovMat();
	}

	private void naiveGWASBed() {
		ChiSquaredDistributionImpl ci = new ChiSquaredDistributionImpl(1);
		SNPFilterPostQC snpPostQC = new SNPFilterPostQC(nssArgs);
		ArrayList<Hukou> hkBook = bed.getHukouBook();

		int[] pIdx = this.data.getMatchedSubjectIdx(1);
		double[][] YY = new double[1+covIdx.length][pIdx.length];


		for (int subjectIdx = 0; subjectIdx < pIdx.length; subjectIdx++) {
			YY[0][subjectIdx] = this.data.getVariable(pheFileIdx, pIdx[subjectIdx], traitIdx[0]);

		}
		ArrayList<NSSGWASResult> nssPheResult = NewIt.newArrayList();
		nssResultList.add(nssPheResult);

		ArrayList<Double> pPheArray = NewIt.newArrayList();
		pArrayList.add(pPheArray);

		for (int i = 0; i < covIdx.length; i++) {
			for (int subjectIdx = 0; subjectIdx < pIdx.length; subjectIdx++) {
				YY[i+1][subjectIdx] = this.data.getVariable(covFileIdx, pIdx[subjectIdx], covIdx[i]);
			}
			ArrayList<NSSGWASResult> nssCovResult = NewIt.newArrayList();
			nssResultList.add(nssCovResult);
			
			ArrayList<Double> pCovArray = NewIt.newArrayList();
			pArrayList.add(pCovArray);
		}

		int numMarkers = map.getMarkerNumberOriginal();
		int numSamples = bed.getNumIndividuals();
		int workingSnpIndex = 0;

		for (int i = 0; i < numMarkers; ++i) {

			if (snpFilter.isSnpIncluded(i)) {
				
				ArrayList<SimpleRegression> sRegArray = NewIt.newArrayList();
				for (int l = 0; l < pArrayList.size(); l++) {
					sRegArray.add(new SimpleRegression());
				}

				int missingCnt = 0;
				int sum = 0;
				int squareSum = 0;

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

						int genotype = (nextByte >> k) & 0b11;

						switch (genotype) {
						case PLINKBinaryParser.HOMOZYGOTE_FIRST:
							break;
						case PLINKBinaryParser.HETEROZYGOTE:
							sum += 1;
							squareSum += 1;
							break;
						case PLINKBinaryParser.HOMOZYGOTE_SECOND:
							sum += 2;
							squareSum += 4;
							break;
						case PLINKBinaryParser.MISSING_GENOTYPE:
							++missingCnt;
							isMissing = true;
							break;
						}
						if (!isMissing) {
							for (int l = 0; l < pArrayList.size(); l++) {
								SimpleRegression sReg = sRegArray.get(l);
								sReg.addData(genotype, YY[l][indCnt - 1]);
							}
						}
					}
				}

				int validSampleCnt = numSamples - missingCnt;
				double variance = 0;
				double average = 0;
				if (validSampleCnt > 2) {
					average = (double) sum / validSampleCnt;
					variance = (squareSum - validSampleCnt * average * average) / (validSampleCnt - 1);
				}

				double freq = 1 - average/2;
				double maf = freq < 0.5 ? freq : (1 - freq);

				boolean isPassPostQC = snpPostQC.isPassPostQC(maf, missingCnt * 1.0D / numSamples);
				if (!isPassPostQC) continue;

				for (int l = 0; l < pArrayList.size(); l++) {
					SimpleRegression sReg = sRegArray.get(l);
					ArrayList<NSSGWASResult> nssResult = nssResultList.get(l);
					ArrayList<Double> pArray = pArrayList.get(l);
					double b, b_se;
					if (freq == 0 || variance == 0) {
						b = Double.NaN;
						b_se = Double.NaN;
					} else {
						b = sReg.getSlope();
						b_se = sReg.getSlopeStdErr();
					}
					SNP snp = map.getSNP(workingSnpIndex);
					NSSGWASResult e1 = new NSSGWASResult(snp, freq, variance, b, b_se);
					nssResult.add(e1);
					pArray.add(e1.GetP());
				}
				workingSnpIndex++;
			}
		}
		snpPostQC.printPostQCSummary();
		
		for (int l = 0; l < pArrayList.size(); l++) {
			ArrayList<Double> pArray = pArrayList.get(l);
			
			Collections.sort(pArray);
			int idx = (int) Math.ceil(pArray.size() / 2);

			Logger.printUserLog("Median of p values is " + pArray.get(idx));

			try {
				double chisq = ci.inverseCumulativeProbability(1 - pArray.get(idx).doubleValue());
				lambdaGC = chisq / 0.4549;
			} catch (MathException e) {
				e.printStackTrace();
			}
			Logger.printUserLog("Lambda GC is: " + (lambdaGC > 0.0001?fmt1.format(lambdaGC):fmt2.format(lambdaGC)));		
			
		}
	}

	private void printCovMat() {
		int[] pIdx = data.getMatchedSubjectIdx(pheFileIdx);
		int[] cIdx = data.getMatchedSubjectIdx(covFileIdx);

		double[][] pheVec = new double[data.getNumberOfSubjects()][traitIdx.length + covIdx.length];

		for (int subjectIdx = 0; subjectIdx < pIdx.length; subjectIdx++) {
			pheVec[subjectIdx][0] = this.data.getVariable(pheFileIdx, pIdx[subjectIdx], this.traitIdx[0]);
			for (int j = 0; j < covIdx.length; j++) {
				pheVec[subjectIdx][j + 1] = this.data.getVariable(covFileIdx, cIdx[subjectIdx], this.covIdx[j]);
			}
		}

		PearsonsCorrelation pc;
		RealMatrix mc;

		pc = new PearsonsCorrelation(pheVec);
		mc = pc.getCorrelationMatrix();

		String Fout = nssArgs.getOutRoot() + ".m.nss";
		PrintStream nssGWAS = FileUtil.CreatePrintStream(Fout);
		for (int i = 0; i < mc.getRowDimension(); i++) {
			for (int j = 0; j < mc.getColumnDimension(); j++) {
				nssGWAS.print(mc.getEntry(i, j) + " ");
			}
			nssGWAS.println();
		}

		nssGWAS.close();
		Logger.printUserLog("");
		Logger.printUserLog("Write the " + mc.getColumnDimension() + "X" + mc.getRowDimension()
				+ " correlation matrix into '" + Fout + "'.");
	}

	public void printResult() {
		
		for (int l = 0; l < pArrayList.size(); l++) {
			String Fout = null;
			if (l == 0) {
				Fout = nssArgs.getOutRoot() + ".p." + (traitIdx[0] + 1) + ".nss";
			} else {
				Fout = nssArgs.getOutRoot() + ".c." + (covIdx[l-1] + 1) + ".nss";
			}
			BufferedWriter nssGZ = FileUtil.ZipFileWriter(Fout+".gz");

			try {
				nssGZ.append(OATHConst.SNP + "\t" + OATHConst.CHR + "\t" + OATHConst.BP + "\t" + OATHConst.RefAle + "\t"
						+ OATHConst.AltAle + "\t" + OATHConst.RAF + "\t" + OATHConst.Vg + "\t" + OATHConst.BETA + "\t"
						+ OATHConst.SE + "\t" + OATHConst.CHI + "\t" + OATHConst.P + "\n");
			} catch (IOException e) {
				Logger.handleException(e,
						"error in writing '" + Fout + ".gz.");
			}

			ArrayList<NSSGWASResult> nssResult = nssResultList.get(l);
			for (int i = 0; i < nssResult.size(); i++) {
				NSSGWASResult e1 = nssResult.get(i);
				try {
					nssGZ.append(e1.printEGWASResult(lambdaGC) +"\n");
				} catch (IOException e) {
					Logger.handleException(e,
							"error in writing '" + Fout + ".gz.");
				}
			}
			try {
				nssGZ.close();
			} catch (IOException e) {
				Logger.handleException(e,
						"error in closing '" + Fout + ".gz.");
			}
			Logger.printUserLog("Write the NSS into '" + Fout + ".gz'.");
		}
	}

	public int getN() {
		return N;
	}
}
