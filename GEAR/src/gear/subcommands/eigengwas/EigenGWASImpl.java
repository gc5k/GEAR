package gear.subcommands.eigengwas;

import gear.data.InputDataSet2;
import gear.family.pedigree.file.MapFile;
import gear.family.pedigree.file.SNP;
import gear.family.plink.PLINKParser;
import gear.family.popstat.GenotypeMatrix;
import gear.family.qc.rowqc.SampleFilter;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.sumstat.qc.rowqc.SumStatQC;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ChiSquaredDistributionImpl;
import org.apache.commons.math.stat.regression.SimpleRegression;

public class EigenGWASImpl extends CommandImpl {
	private EigenGWASArguments eigenArgs;
	private MapFile mapFile;
	private SampleFilter sf;
	private SumStatQC ssQC;
	private GenotypeMatrix gm;
	private int traitIdx;
	private InputDataSet2 data = new InputDataSet2();
	private ArrayList<EigenGWASResult> eGWASResult = NewIt.newArrayList();

	private double lambdaGC = 1;
	private int monoLoci = 0;

	public void execute(CommandArguments cmdArgs) 
	{
		this.eigenArgs = (EigenGWASArguments) cmdArgs;

		PLINKParser pp = PLINKParser.parse(this.eigenArgs);
		this.sf = new SampleFilter(pp.getPedigreeData(), pp.getMapData());
		this.ssQC = new SumStatQC(pp.getPedigreeData(), pp.getMapData(), this.sf);
		this.mapFile = this.ssQC.getMapFile();
		this.gm = new GenotypeMatrix(this.ssQC.getSample());
		
		this.traitIdx = this.eigenArgs.getMpheno()[0];
		this.data.addFile(this.eigenArgs.getFam());
		this.data.addFile(this.eigenArgs.getPhenotypeFile(), this.eigenArgs.getMpheno());
		if (this.eigenArgs.getKeepFile() != null)
		{
			this.data.addFile(this.eigenArgs.getKeepFile());
		}
		data.LineUpFiles();


		eigenGWAS();
		printResult();
	}

	private void eigenGWAS()
	{
		ChiSquaredDistributionImpl ci = new ChiSquaredDistributionImpl(1);
		ArrayList<SNP> snpList = this.mapFile.getMarkerList();
		
		int[] gIdx = this.data.getMatchedSubjectIdx(0);
		int[] pIdx = this.data.getMatchedSubjectIdx(1);
		
		double[] Y = new double[pIdx.length];
		ArrayList<Double> pArray = NewIt.newArrayList();
		double threshold = 0.0D;

		for (int subjectIdx = 0; subjectIdx < Y.length; subjectIdx++) 
		{
			Y[subjectIdx] = this.data.getVariable(1, pIdx[subjectIdx], this.traitIdx);
			threshold += Y[subjectIdx];
		}
		threshold /= Y.length;

		for (int i = 0; i < this.gm.getNumMarker(); i++)
		{
			SNP snp = (SNP) snpList.get(i);

			if ((!this.eigenArgs.isChrFlagOn()) || (Integer.parseInt(snp.getChromosome()) == this.eigenArgs.getChr())) 
			{
				SimpleRegression sReg = new SimpleRegression();
				double n1 = 0.0D;
				double n2 = 0.0D;
				double N = 0.0D;
				double freq1 = 0.0D;
				double freq2 = 0.0D;
				double freq = 0.0D;

				for (int j = 0; j < pIdx.length; j++)
				{
					int g = this.gm.getAdditiveScoreOnFirstAllele(gIdx[j], i);
					if (g != 3) {
						sReg.addData(g, Y[j]);
						if (Y[j] < threshold) 
						{
							n1 += 1.0D;
							freq1 += g;
						}
						else
						{
							n2 += 1.0D;
							freq2 += g;
						}
						N += 1.0D;
						freq += g;
					}
				}
				if(snp.isMonopolic() && N <=1) 
				{
					monoLoci++;
					continue;
				}
				freq1 /= 2.0D * n1;
				freq2 /= 2.0D * n2;
				freq /= 2.0D * N;

				double fst = 2 * (n1 / N * (freq1 - freq) * (freq1 - freq) + n2 / N * (freq2 - freq) * (freq2 - freq))
						/ (freq * (1.0D - freq));

				double b = sReg.getSlope();
				double b_se = sReg.getSlopeStdErr();
				EigenGWASResult e1 = new EigenGWASResult(snp, freq, b, b_se, n1, freq1, n2, freq2, fst);
				eGWASResult.add(e1);
				pArray.add(e1.GetP());
			}
		}
		Collections.sort(pArray);
		int idx = (int) Math.ceil(pArray.size() / 2);

		if (monoLoci > 1)
		{
			Logger.printUserLog("Removed " + monoLoci + " monomorphic loci.");			
		}
		else if (monoLoci == 1)
		{
			Logger.printUserLog("Removed " + monoLoci + " monomorphic locus.");
		}

		Logger.printUserLog("Median of p values is " + pArray.get(idx));

		try 
		{
			double chisq = ci.inverseCumulativeProbability(1 - pArray.get(idx).doubleValue());
			lambdaGC = chisq / 0.4549;
		} 
		catch (MathException e) 
		{
			e.printStackTrace();
		}

		Logger.printUserLog("Lambda GC is : " + lambdaGC);
	}

	public void printResult() 
	{

		PrintStream eGWAS = FileUtil.CreatePrintStream(this.eigenArgs.getOutRoot() + ".egwas");
		eGWAS.println("SNP\tCHR\tBP\tRefAllele\tAltAllele\tfreq\tBeta\tSE\tChi\tP\tPGC\tn1\tfreq1\tn2\tfreq2\tFst");

		for (int i = 0; i < eGWASResult.size(); i++) 
		{
			EigenGWASResult e1 = eGWASResult.get(i);
			eGWAS.println(e1.printEGWASResult(lambdaGC));
		}
		eGWAS.close();
	}
}