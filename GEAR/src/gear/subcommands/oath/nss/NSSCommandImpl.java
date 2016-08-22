package gear.subcommands.oath.nss;

import gear.data.InputDataSet;
import gear.family.pedigree.file.MapFile;
import gear.family.pedigree.file.SNP;
import gear.family.plink.PLINKParser;
import gear.family.popstat.GenotypeMatrix;
import gear.family.qc.rowqc.SampleFilter;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.subcommands.oath.OATHConst;
import gear.sumstat.qc.rowqc.SumStatQC;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ChiSquaredDistributionImpl;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.stat.regression.SimpleRegression;
import org.apache.commons.math.stat.StatUtils;
import org.apache.commons.math.stat.correlation.PearsonsCorrelation;

public class NSSCommandImpl extends CommandImpl 
{
	private NSSCommandArguments nssArgs;
	private MapFile mapFile;
	private SampleFilter sf;
	private SumStatQC ssQC;
	private GenotypeMatrix gm;
	private int[] traitIdx;
	private InputDataSet data = null;
	private ArrayList<NSSGWASResult> nssResult;

	private double lambdaGC = 1;
	private int monoLoci = 0;

	public void execute(CommandArguments cmdArgs) 
	{
		this.nssArgs = ((NSSCommandArguments) cmdArgs);

		this.traitIdx = this.nssArgs.getMpheno();
		int[] covIdx = new int[this.traitIdx.length - 1];
		System.arraycopy(this.traitIdx, 1, covIdx, 0, covIdx.length);
		this.data = new InputDataSet(this.nssArgs.getFam(), this.nssArgs.getPhenotypeFile(), this.nssArgs.getPhenotypeFile(), this.nssArgs.getMpheno()[0], covIdx);
//		this.data.readSubjectIDFile(this.nssArgs.getFam());
//		this.data.readPhenotypeFile(this.nssArgs.getPhenotypeFile());

		PLINKParser pp = PLINKParser.parse(this.nssArgs);
		this.sf = new SampleFilter(pp.getPedigreeData(), pp.getMapData());
		this.ssQC = new SumStatQC(pp.getPedigreeData(), pp.getMapData(), this.sf);
		this.mapFile = this.ssQC.getMapFile();
		this.gm = new GenotypeMatrix(this.ssQC.getSample());

		for (int i = 0; i < traitIdx.length; i++)
		{
			Logger.printUserLog("");
			Logger.printUserLog("Generating naive summary statistics (NSS) for "+ (traitIdx[i] + 1) + "th variable.");
			nssResult = NewIt.newArrayList();
			naiveGWAS(i);
			printResult(i);
		}

		String Fout = nssArgs.getOutRoot()  + ".list.nss";
		PrintStream nssList = FileUtil.CreatePrintStream(Fout);
		for (int i = 0; i < traitIdx.length; i++)
		{
			nssList.println(nssArgs.getOutRoot() + "." + (traitIdx[i] + 1) + ".nss");
		}
		nssList.close();
		printCovMat();
	}

	private void printCovMat()
	{
		ArrayList<ArrayList<Double>> Dat = NewIt.newArrayList();
		int[] pheIdx = data.getMatchedPheSubIdx();
		for (int subjectIdx = 0; subjectIdx < pheIdx.length; subjectIdx++)
		{
			ArrayList<Double> dat = NewIt.newArrayList();
			boolean isMissing = false;
			for (int j = 0; j < traitIdx.length; j++)
			{
				if (!this.data.isPhenotypeMissing(pheIdx[subjectIdx], this.traitIdx[j]))
				{
					dat.add((double) this.data.getPhenotype(pheIdx[subjectIdx], this.traitIdx[j]));
				}
				else
				{
					isMissing = true;
				}
			}
			if (!isMissing)
			{
				Dat.add(dat);
			}
		}

		double[][] pheVec = new double[Dat.size()][traitIdx.length];
		for (int i = 0; i < Dat.size(); i++)
		{
			ArrayList<Double> D = Dat.get(i);
			for (int j = 0; j < D.size(); j++)
			{
				pheVec[i][j] = D.get(j);
			}
		}

		PearsonsCorrelation pc = new PearsonsCorrelation(pheVec);
		RealMatrix mc = pc.getCorrelationMatrix();

		String Fout = nssArgs.getOutRoot()  + ".m.nss";
		PrintStream nssGWAS = FileUtil.CreatePrintStream(Fout);
		for(int i = 0; i < mc.getRowDimension(); i++)
		{
			for(int j = 0; j < mc.getColumnDimension(); j++)
			{
				nssGWAS.print(mc.getEntry(i, j) + " ");
			}
			nssGWAS.println();
		}

		nssGWAS.close();
		Logger.printUserLog("");
		Logger.printUserLog("Write the " + mc.getColumnDimension() + "X" + mc.getRowDimension() + " correlation matrix into '" + Fout + "'.");
	}

	private void naiveGWAS(int tIdx) 
	{
		ChiSquaredDistributionImpl ci = new ChiSquaredDistributionImpl(1);
		ArrayList<SNP> snpList = this.mapFile.getMarkerList();

		int[] pIdx = data.getMatchedPheSubIdx();
		double[] Y = new double[pIdx.length];
		ArrayList<Double> pArray = NewIt.newArrayList();

		for (int subjectIdx = 0; subjectIdx < Y.length; subjectIdx++) 
		{
			Y[subjectIdx] = this.data.getPhenotype(pIdx[subjectIdx], this.traitIdx[tIdx]);
		}
		Y = StatUtils.normalize(Y);

		int[] subIdx = data.getMatchedSubIdx();
		for (int i = 0; i < this.gm.getNumMarker(); i++)
		{
			SNP snp = (SNP) snpList.get(i);

			if ((!this.nssArgs.isChrFlagOn()) || (Integer.parseInt(snp.getChromosome()) == this.nssArgs.getChr()))
			{
				SimpleRegression sReg = new SimpleRegression();
				double N = 0.0D;
				double freq = 0.0D;
				double xx = 0;
				double mx = 0;

				for (int j = 0; j < subIdx.length; j++)
				{
					int idx = subIdx[j];
					int g = this.gm.getAdditiveScoreOnFirstAllele(idx, i);
					if (g != 3) {
						sReg.addData(g, Y[idx]);
						N += 1.0D;
						freq += g;
						mx += g;
						xx += g*g;
					}
				}

				if (snp.isMonopolic() && N <= 1) 
				{
					monoLoci++;
					continue;
				}

				freq /= 2.0D * N;
				double vg = (xx - (mx/N) * (mx/N) * N)/(N-1);

				double b = sReg.getSlope();
				double b_se = sReg.getSlopeStdErr();

				NSSGWASResult e1 = new NSSGWASResult(snp, freq, vg, b, b_se);
				nssResult.add(e1);
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
		Logger.printUserLog("Lambda GC is: " + lambdaGC);
	}

	public void printResult(int tIdx)
	{
		int fidx = this.traitIdx[tIdx] + 1;
		String Fout = nssArgs.getOutRoot()  + "." + fidx + ".nss";
		PrintStream nssGWAS = FileUtil.CreatePrintStream(Fout);
		nssGWAS.println(OATHConst.SNP +"\t" + OATHConst.CHR + "\t" + OATHConst.BP+ "\t" + OATHConst.RefAle + "\t" + OATHConst.AltAle +"\t" + OATHConst.Freq + "\t" + OATHConst.Vg + "\t" + OATHConst.BETA +"\t" + OATHConst.SE + "\t" + OATHConst.CHI + "\t" + OATHConst.P);

		for (int i = 0; i < nssResult.size(); i++)
		{
			NSSGWASResult e1 = nssResult.get(i);
			nssGWAS.println(e1.printEGWASResult(lambdaGC));
		}
		nssGWAS.close();
		
		Logger.printUserLog("Write the NSS into '" + Fout + "'.");
	}
}
