package gear.subcommands.eigengwasepi;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ChiSquaredDistributionImpl;
import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;

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
import gear.util.pop.PopStat;

public class EigenGWASEpiCommandImpl extends CommandImpl
{

	private EigenGWASEpiCommandArguments eigenArgs;
	private MapFile mapFile;
	private SampleFilter sf;
	private SumStatQC ssQC;
	private GenotypeMatrix gm;
	private double[][] gfreq;
	private int traitIdx;
	private InputDataSet2 data = new InputDataSet2();
	private ArrayList<EigenGWASEpiResult> EpiGWASResult = NewIt.newArrayList();

	private double lambdaGC_A1 = 1;
	private double lambdaGC_D1 = 1;
	private double lambdaGC_A2 = 1;
	private double lambdaGC_D2 = 1;
	private double lambdaGC_AA = 1;

	private int monoLoci = 0;
	private int singularLoci = 0;

	public void execute(CommandArguments cmdArgs) 
	{
		this.eigenArgs = (EigenGWASEpiCommandArguments) cmdArgs;

		this.traitIdx = this.eigenArgs.getMpheno()[0];
		data.addFile(this.eigenArgs.getFam());
		data.addFile(this.eigenArgs.getPhenotypeFile(), this.eigenArgs.getMpheno());
		if (eigenArgs.getKeepFile() != null)
		{
			data.addFile(this.eigenArgs.getKeepFile());
		}
		data.LineUpFiles();

		PLINKParser pp = PLINKParser.parse(this.eigenArgs);
		this.sf = new SampleFilter(pp.getPedigreeData(), pp.getMapData());
		this.ssQC = new SumStatQC(pp.getPedigreeData(), pp.getMapData(), this.sf);
		this.mapFile = this.ssQC.getMapFile();
		this.gm = new GenotypeMatrix(this.ssQC.getSample());

		eigenEpiGWAS();
		printResult();
	}

	private void eigenEpiGWAS()
	{
		ChiSquaredDistributionImpl ci = new ChiSquaredDistributionImpl(1);
		ArrayList<SNP> snpList = this.mapFile.getMarkerList();

		int[] gIdx = this.data.getMatchedSubjectIdx(0);
		int[] pIdx = this.data.getMatchedSubjectIdx(1);

		double[] Y = new double[pIdx.length];

		ArrayList<Double> pA1Array = NewIt.newArrayList();
		ArrayList<Double> pD1Array = NewIt.newArrayList();
		ArrayList<Double> pA2Array = NewIt.newArrayList();
		ArrayList<Double> pD2Array = NewIt.newArrayList();
		ArrayList<Double> pAAArray = NewIt.newArrayList();

		double threshold = 0.0D;

		for (int subjectIdx = 0; subjectIdx < Y.length; subjectIdx++) 
		{
			Y[subjectIdx] = this.data.getVariable(1, pIdx[subjectIdx], this.traitIdx);
			threshold += Y[subjectIdx];
		}
		threshold /= Y.length;

		// PrintStream eGWAS =
		// FileUtil.CreatePrintStream(this.eigenArgs.getOutRoot() + ".egwas");
		// eGWAS.println("SNP\tCHR\tBP\tRefAllele\tAltAllele\tfreq\tBeta\tSE\tChi\tP\tPGC\tn1\tfreq1\tn2\tfreq2\tFst");

		gfreq = PopStat.calGenoFrequency(gm, gm.getNumMarker());
		for (int i = 0; i < this.gm.getNumMarker(); i++)
		{
			SNP snp1 = (SNP) snpList.get(i);

			double[][] x= new double[gIdx.length][4]; //a1+d1+a2+d2
				if (snp1.isMonopolic())
				{
					monoLoci++;
					continue;
				}

				double N1_1 = 0.0D;
				double N1_2 = 0.0D;
				double N1 = 0.0D;
				double Freq1_1 = 0.0D;
				double Freq1_2 = 0.0D;
				double Freq1 = 0.0D;
				for (int j = 0; j < gIdx.length; j++)
				{
					double[] cd1 = getCoding(gIdx[j], i);
					if (cd1[2] == 1)
					{
						int g = this.gm.getAdditiveScoreOnFirstAllele(gIdx[j], i);
						x[j][0] = cd1[0];
						x[j][1] = cd1[1];

						if (Y[j] < threshold)
						{
							N1_1 += 1.0D;
							Freq1_1 += g;
						}
						else
						{
							N1_2 += 1.0D;
							Freq1_2 += g;
						}
						N1 += 1.0D;
						Freq1 += g;
					}
				}

				Freq1_1 /= 2.0D * N1_1;
				Freq1_2 /= 2.0D * N1_2;
				Freq1 /= 2.0D * N1;

				double fst1 = 2 * (N1_1 / N1 * (Freq1_1 - Freq1) * (Freq1_1 - Freq1) + N1_2 / N1 * (Freq1_2 - Freq1) * (Freq1_2 - Freq1))
						/ (Freq1 * (1.0D - Freq1));

				for (int i2 = i+1; i2 < this.gm.getNumMarker(); i2++)
				{
					SNP snp2 = (SNP) snpList.get(i2);

					if (snp2.isMonopolic())
					{
						monoLoci++;
						continue;
					}

					double N2_1 = 0.0D;
					double N2_2 = 0.0D;
					double N2 = 0.0D;
					double Freq2_1 = 0.0D;
					double Freq2_2 = 0.0D;
					double Freq2 = 0.0D;

					for (int j2 = 0; j2 < gIdx.length; j2++)
					{
						double[] cd2 = getCoding(gIdx[j2], i2);

						if (cd2[2] == 1)
						{
							int g2 = this.gm.getAdditiveScoreOnFirstAllele(gIdx[j2], i);

							x[j2][2] = cd2[0];
							x[j2][3] = cd2[1];

							if (Y[j2] < threshold)
							{
								N2_1 += 1.0D;
								Freq2_1 += g2;
							}
							else
							{
								N2_2 += 1.0D;
								Freq2_2 += g2;
							}
							N2 += 1.0D;
							Freq2 += g2;
						}
					}

					Freq2_1 /= 2.0D * N2_1;
					Freq2_2 /= 2.0D * N2_2;
					Freq2 /= 2.0D * N2;

					double fst2 = 2 * (N2_1 / N2 * (Freq2_1 - Freq2) * (Freq2_1 - Freq2) + N2_2 / N2 * (Freq2_2 - Freq2) * (Freq2_2 - Freq2))
							/ (Freq2 * (1.0D - Freq2));

					double[][] x1= new double[gIdx.length][6];//u + a1 + d1 + a2 + d2 + aa
					for (int k1 = 0; k1 < x1.length; k1++)
					{
						x1[k1][0] = 1;
						x1[k1][1] = x[k1][0];
						x1[k1][2] = x[k1][1];
						x1[k1][3] = x[k1][2];
						x1[k1][4] = x[k1][3];
						x1[k1][5] = x[k1][0]*x[k1][2];
					}

					RealMatrix g1 = new Array2DRowRealMatrix(x1);

					if (eigenArgs.isInbred())
					{
						int[] gidx = new int[g1.getRowDimension()];
						for(int g_i = 0; g_i < g1.getRowDimension(); g_i++) gidx[g_i] = g_i;
						int[] pidx = {0,1,3,5};
						g1 = g1.getSubMatrix(gidx, pidx);
					}

					RealMatrix XtX = g1.transpose().multiply(g1);

					boolean isNonSingular = (new LUDecompositionImpl(XtX)).getSolver().isNonSingular();

					if (!isNonSingular)
					{
						Logger.printUserLog("Model is singular for '" + snp1.getName() + "' + '" + snp2.getName() + "'; "+ "skipped.");
						singularLoci++;
						continue;
					}
					else
					{
						RealMatrix XtX_inv = (new LUDecompositionImpl(XtX)).getSolver().getInverse();

						RealMatrix g1_tran=g1.transpose();
						RealMatrix y = new Array2DRowRealMatrix(Y);
						RealMatrix B = XtX_inv.multiply(g1_tran).multiply(y);
						RealMatrix SST = y.transpose().multiply(y);
						RealMatrix SSR = B.transpose().multiply(g1_tran).multiply(y);
						double sse = (SST.getEntry(0, 0) - SSR.getEntry(0, 0))/(y.getRowDimension() - B.getRowDimension());

						RealMatrix BV = XtX_inv.scalarMultiply(sse);
						if (!eigenArgs.isInbred())
						{
							if(BV.getEntry(1, 1) > 0 && BV.getEntry(2, 2) > 0 && BV.getEntry(3, 3) > 0 && BV.getEntry(4, 4) > 0 && BV.getEntry(5, 5) > 0)
							{
								EigenGWASEpiResult e1 = new EigenGWASEpiResult(snp1, Freq1, fst1, B.getEntry(1, 0), Math.sqrt(BV.getEntry(1, 1)), B.getEntry(2, 0), Math.sqrt(BV.getEntry(2, 2)), 
										snp2, Freq2, fst2, B.getEntry(3, 0), Math.sqrt(BV.getEntry(3, 3)), B.getEntry(4, 0), Math.sqrt(BV.getEntry(4, 4)), B.getEntry(5, 0), Math.sqrt(BV.getEntry(5, 5)));
								EpiGWASResult.add(e1);

								pA1Array.add(e1.GetAP1());
								pD1Array.add(e1.GetDP1());
								pA2Array.add(e1.GetAP2());
								pD2Array.add(e1.GetDP2());
								pAAArray.add(e1.GetAAP());
							}							
						}
						else
						{
							if(BV.getEntry(1, 1) > 0 && BV.getEntry(2, 2) > 0 && BV.getEntry(3, 3) > 0)
							{
								EigenGWASEpiResult e1 = new EigenGWASEpiResult(snp1, Freq1, fst1, B.getEntry(1, 0), Math.sqrt(BV.getEntry(1, 1)), 1, 1, 
										snp2, Freq2, fst2, B.getEntry(2, 0), Math.sqrt(BV.getEntry(2, 2)), 1, 1, B.getEntry(3, 0), Math.sqrt(BV.getEntry(3, 3)));
								EpiGWASResult.add(e1);

								pA1Array.add(e1.GetAP1());
								pD1Array.add(e1.GetDP1());
								pA2Array.add(e1.GetAP2());
								pD2Array.add(e1.GetDP2());
								pAAArray.add(e1.GetAAP());
							}
						}
				}
			}
		}

		Collections.sort(pA1Array);
		int idxA1 = (int) Math.ceil(pA1Array.size() / 2);
		Collections.sort(pD1Array);
		int idxD1 = (int) Math.ceil(pD1Array.size() / 2);
		Collections.sort(pA2Array);
		int idxA2 = (int) Math.ceil(pA2Array.size() / 2);
		Collections.sort(pD2Array);
		int idxD2 = (int) Math.ceil(pD2Array.size() / 2);
		Collections.sort(pAAArray);
		int idxAA = (int) Math.ceil(pAAArray.size() / 2);

		if (monoLoci > 1)
		{
			Logger.printUserLog("Removed " + monoLoci + " monomorphic loci.");			
		}
		else if (monoLoci == 1)
		{
			Logger.printUserLog("Removed " + monoLoci + " monomorphic locus.");
		}

		if (singularLoci > 1)
		{
			Logger.printUserLog("Removed " + singularLoci + " singular loci.");
		}
		else if (singularLoci == 1)
		{
			Logger.printUserLog("Removed " + singularLoci + " singular locus.");			
		}
//		Logger.printUserLog("Median of p values (additive 1) is " + pA1Array.get(idxA1));
//		Logger.printUserLog("Median of p values (dominance 1) is " + pD1Array.get(idxD1));
//		Logger.printUserLog("Median of p values (additive 2) is " + pA2Array.get(idxA2));
//		Logger.printUserLog("Median of p values (dominance 2) is " + pD2Array.get(idxD2));
		Logger.printUserLog("Median of p values (AA epi) is " + pD2Array.get(idxD2));

		try
		{
			double chisqA1 = ci.inverseCumulativeProbability(1 - pA1Array.get(idxA1).doubleValue());
			lambdaGC_A1 = chisqA1 / 0.4549;

			double chisqD1 = ci.inverseCumulativeProbability(1 - pD1Array.get(idxD1).doubleValue());
			lambdaGC_D1 = chisqD1 / 0.4549;
			
			double chisqA2 = ci.inverseCumulativeProbability(1 - pA2Array.get(idxA2).doubleValue());
			lambdaGC_A2 = chisqA2 / 0.4549;

			double chisqD2 = ci.inverseCumulativeProbability(1 - pD2Array.get(idxD2).doubleValue());
			lambdaGC_D2 = chisqD2 / 0.4549;

			double chisqAA = ci.inverseCumulativeProbability(1 - pAAArray.get(idxAA).doubleValue());
			lambdaGC_AA = chisqAA / 0.4549;
		}
		catch (MathException e)
		{
			e.printStackTrace();
		}

//		Logger.printUserLog("Lambda GC (additive 1) is: " + lambdaGC_A1);
//		Logger.printUserLog("Lambda GC (dominance 1) is: " + lambdaGC_D1);
//		Logger.printUserLog("Lambda GC (additive 2) is: " + lambdaGC_A2);
//		Logger.printUserLog("Lambda GC (dominance 2) is: " + lambdaGC_D2);
		Logger.printUserLog("Lambda GC (AA) is: " + lambdaGC_AA);		
	}

	private double[] getCoding(int ind, int locus)
	{
		double[] cd = new double[3];
		int g = this.gm.getAdditiveScoreOnFirstAllele(ind, locus);
		if (g != 3)
		{
			cd[2] = 1;
			if (g == 0)
			{
				/*
				x[j][0] = -1 * (2 * gfreq[i][2] + gfreq[i][1]);
				x[j][1] = 1/ (8 * gfreq[i][0]);
				*/
				cd[0] = -gfreq[locus][1] - 2*gfreq[locus][2];
				cd[1] = (2 * gfreq[locus][1] * gfreq[locus][2])/(gfreq[locus][0]+gfreq[locus][2]-(gfreq[locus][0]-gfreq[locus][2])*(gfreq[locus][0]-gfreq[locus][2]));
			}
			else if (g == 1)
			{
				/*
				x[j][0] = gfreq[i][0] - gfreq[i][2];
				x[j][1] = -1 / (4 * gfreq[i][1]);
				*/
				cd[0] = 1 - gfreq[locus][1] - 2*gfreq[locus][2];
				cd[1] = -1 * (4*gfreq[locus][0]*gfreq[locus][2])/(gfreq[locus][0]+gfreq[locus][2]-(gfreq[locus][0]-gfreq[locus][2])*(gfreq[locus][0]-gfreq[locus][2]));
			}
			else
			{
				/*
				x[j][0] = -1 * (2 - 1*gfreq[i][1] - 2*gfreq[i][2]);
				x[j][1] = -1 / (8 * gfreq[i][2]);
				*/
				cd[0] = 2 - gfreq[locus][1] - 2*gfreq[locus][2];
				cd[1] = (2*gfreq[locus][0]*gfreq[locus][1])/(gfreq[locus][0]+gfreq[locus][2]-(gfreq[locus][0]-gfreq[locus][2])*(gfreq[locus][0]-gfreq[locus][2]));
			}
		}
		return cd;
	}
	
	public void printResult() 
	{

		PrintStream eGWAS = FileUtil.CreatePrintStream(this.eigenArgs.getOutRoot() + ".egwasepi");
		eGWAS.println("SNP1\tCHR1\tBP1\tRefAllele1\tAltAllele1\tFreq1\tBeta1\tSE1\tChi1\tP1\tDom1\tDomSE1\tChisqDom1\tpDom1\tFst1\tSNP2\tCHR2\tBP2\tRefAllele2\tAltAllele2\tFreq2\tBeta2\tSE2\tChi2\tP2\tDom2\tDomSE2\tChisqDom2\tpDom2\tFst2\tAA\tAASE\tChisqAA\tpAA");

		for (int i = 0; i < EpiGWASResult.size(); i++)
		{
			EigenGWASEpiResult e1 = EpiGWASResult.get(i);
			eGWAS.println(e1.printEGWASEpiResult());
		}
		eGWAS.close();
	}

}
