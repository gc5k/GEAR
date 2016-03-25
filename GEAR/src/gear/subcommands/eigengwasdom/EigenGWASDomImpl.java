package gear.subcommands.eigengwasdom;

import gear.data.InputDataSet;
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
import gear.util.pop.PopStat;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ChiSquaredDistributionImpl;
import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;

public class EigenGWASDomImpl extends CommandImpl 
{
	private EigenGWASDomCommandArguments eigenArgs;
	private MapFile mapFile;
	private SampleFilter sf;
	private SumStatQC ssQC;
	private GenotypeMatrix gm;
	private int traitIdx;
	private InputDataSet data = new InputDataSet();
	private ArrayList<EigenGWASDomResult> eGWASResult = NewIt.newArrayList();

	private double lambdaGC = 1;
	private int monoLoci = 0;
	private int singularLoci = 0;

	public void execute(CommandArguments cmdArgs) 
	{
		this.eigenArgs = ((EigenGWASDomCommandArguments) cmdArgs);

		this.traitIdx = this.eigenArgs.getMpheno();
		this.data.readSubjectIDFile(this.eigenArgs.getFam());
		this.data.readPhenotypeFile(this.eigenArgs.getPhenotypeFile());

		PLINKParser pp = PLINKParser.parse(this.eigenArgs);
		this.sf = new SampleFilter(pp.getPedigreeData(), pp.getMapData());
		this.ssQC = new SumStatQC(pp.getPedigreeData(), pp.getMapData(), this.sf);
		this.mapFile = this.ssQC.getMapFile();
		this.gm = new GenotypeMatrix(this.ssQC.getSample());

		eigenGWAS();
		printResult();
	}

	private void eigenGWAS() 
	{
		ChiSquaredDistributionImpl ci = new ChiSquaredDistributionImpl(1);
		ArrayList<SNP> snpList = this.mapFile.getMarkerList();

		double[] Y = new double[this.data.getNumberOfSubjects()];
		ArrayList<Integer> pheIdx = NewIt.newArrayList();
		ArrayList<Double> pArray = NewIt.newArrayList();
		double threshold = 0.0D;

		for (int subjectIdx = 0; subjectIdx < Y.length; subjectIdx++) 
		{
			if (!this.data.isPhenotypeMissing(subjectIdx, this.traitIdx))
			{
				pheIdx.add(Integer.valueOf(subjectIdx));
				Y[subjectIdx] = this.data.getPhenotype(subjectIdx, this.traitIdx);
				threshold += Y[subjectIdx];
			}
		}
		threshold /= pheIdx.size();

		// PrintStream eGWAS =
		// FileUtil.CreatePrintStream(this.eigenArgs.getOutRoot() + ".egwas");
		// eGWAS.println("SNP\tCHR\tBP\tRefAllele\tAltAllele\tfreq\tBeta\tSE\tChi\tP\tPGC\tn1\tfreq1\tn2\tfreq2\tFst");

		double[][] gfreq = PopStat.calGenoFrequency(gm, gm.getNumMarker());
		for (int i = 0; i < this.gm.getNumMarker(); i++)
		{
			SNP snp = (SNP) snpList.get(i);
			if (snp.isMonopolic())
			{
				monoLoci++;
				continue;
			}

			double[][] x= new double[this.gm.getNumIndivdial()][2];
			if ((!this.eigenArgs.isChrFlagOn()) || (Integer.parseInt(snp.getChromosome()) == this.eigenArgs.getChr())) 
			{
				double n1 = 0.0D;
				double n2 = 0.0D;
				double N = 0.0D;
				double freq1 = 0.0D;
				double freq2 = 0.0D;
				double freq = 0.0D;
				for (int j = 0; j < pheIdx.size(); j++)
				{
					int idx = ((Integer) pheIdx.get(j)).intValue();
					int g = this.gm.getAdditiveScoreOnFirstAllele(idx, i);
					if (g != 3)
					{
						if (g == 0)
						{
							x[j][0] = -1 * (2 * gfreq[i][2] + gfreq[i][1]);
							x[j][1] = 1/ (8 * gfreq[i][0]);
						}
						else if (g == 1)
						{
							x[j][0] = gfreq[i][0] - gfreq[i][2];
							x[j][1] = -1 / (4 * gfreq[i][1]);	
						}
						else
						{
							x[j][0] = (2 * gfreq[i][0] + gfreq[i][1]);
							x[j][1] = -1 / (8 * gfreq[i][2]);			
						}

						if (Y[idx] < threshold)
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

				freq1 /= 2.0D * n1;
				freq2 /= 2.0D * n2;
				freq /= 2.0D * N;

				double fst = 2 * (n1 / N * (freq1 - freq) * (freq1 - freq) + n2 / N * (freq2 - freq) * (freq2 - freq))
						/ (freq * (1.0D - freq));

				double[][] x1= new double[this.gm.getNumIndivdial()][3];
				for(int k1 = 0; k1 < x1.length; k1++)
				{
					x1[k1][0] = 1;
					x1[k1][1] = x[k1][0];
					x1[k1][2] = x[k1][1];
				}

				RealMatrix g1 = new Array2DRowRealMatrix(x1);
				RealMatrix XtX = g1.transpose().multiply(g1);

				boolean isNonSingular = (new LUDecompositionImpl(XtX)).getSolver().isNonSingular();

				if (!isNonSingular)
				{
					Logger.printUserLog("Model is singular for " + snp.getName() + ".");
					singularLoci++;
					continue;
				}
				else
				{
					System.out.println(XtX);
					RealMatrix XtX_inv = (new LUDecompositionImpl(XtX)).getSolver().getInverse();
					System.out.println(XtX_inv);

					RealMatrix g1_tran=g1.transpose();
					RealMatrix y = new Array2DRowRealMatrix(Y);
					RealMatrix B = XtX_inv.multiply(g1_tran).multiply(y);
					System.out.println(B);
					RealMatrix SST = y.transpose().multiply(y);
					System.out.println(SST);
					RealMatrix SSR = B.transpose().multiply(g1_tran).multiply(y);
					System.out.println(SSR);
					double sse = (SST.getEntry(0, 0) - SSR.getEntry(0, 0))/(y.getRowDimension() - B.getRowDimension());

					RealMatrix BV = XtX_inv.scalarMultiply(sse);
					System.out.println(BV);
					if(BV.getEntry(1, 1) > 0 && BV.getEntry(2, 2) > 0)
					{
						EigenGWASDomResult e1 = new EigenGWASDomResult(snp, freq, B.getEntry(1, 0), Math.sqrt(BV.getEntry(1, 1)), B.getEntry(2, 0), Math.sqrt(BV.getEntry(2, 2)), n1, freq1, n2, freq2, fst);
						eGWASResult.add(e1);
						pArray.add(e1.GetP());
					}
				}

//				mReg.newSampleData(Y, x);
//				double[] beta = mReg.estimateRegressionParameters();
//				double[] bV = mReg.estimateRegressionParametersStandardErrors();

//				System.out.println(fst + " " + gfreq[i][0] + " " + gfreq[i][1] + " " + gfreq[i][2]);
//				for(int k = 0; k < x.length; k++)
//				{
//					System.out.println(x[k][0] + " " + x[k][1]);
//				}
//

//				System.exit(1);

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
		eGWAS.println("SNP\tCHR\tBP\tRefAllele\tAltAllele\tfreq\tBeta\tSE\tChi\tP\tPGC\tDom\tDomSE\tpDom\tn1\tfreq1\tn2\tfreq2\tFst");

		for (int i = 0; i < eGWASResult.size(); i++) 
		{
			EigenGWASDomResult e1 = eGWASResult.get(i);
			eGWAS.println(e1.printEGWASResult(lambdaGC));
		}
		eGWAS.close();
	}
	
}