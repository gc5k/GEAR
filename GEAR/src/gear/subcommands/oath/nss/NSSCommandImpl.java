package gear.subcommands.oath.nss;

import gear.data.InputDataSet2;
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
import gear.util.pop.PopStat;

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
	private int[] covIdx;
	private double[][] frq;
	private int N;
	private InputDataSet2 data = null;
	private ArrayList<NSSGWASResult> nssResult;

	private double lambdaGC = 1;
//	private int monoLoci = 0;
	
	private int famFileIdx = 0;
	private int pheFileIdx = 1;
	private int covFileIdx = 2;

	public void execute(CommandArguments cmdArgs) 
	{
		this.nssArgs = ((NSSCommandArguments) cmdArgs);

		this.traitIdx = this.nssArgs.getMpheno();
		this.covIdx = this.nssArgs.getCovNumber();

		data = new InputDataSet2();
		data.addFile(this.nssArgs.getFam());
		data.addFile(this.nssArgs.getPhenotypeFile(), this.nssArgs.getMpheno());
		data.addFile(this.nssArgs.getCovFile(), this.nssArgs.getCovNumber());
		if(this.nssArgs.getKeepFile() != null) data.addFile(this.nssArgs.getKeepFile());
		data.LineUpFiles();

		this.N = data.getNumberOfSubjects();

		PLINKParser pp = PLINKParser.parse(this.nssArgs);
		this.sf = new SampleFilter(pp.getPedigreeData(), pp.getMapData());
		this.ssQC = new SumStatQC(pp.getPedigreeData(), pp.getMapData(), this.sf);
		this.mapFile = this.ssQC.getMapFile();
		this.gm = new GenotypeMatrix(this.ssQC.getSample());

		frq = PopStat.calAlleleFrequency(this.gm, this.gm.getNumMarker());
		for (int i = 0; i < traitIdx.length; i++)
		{
			Logger.printUserLog("");
			Logger.printUserLog("Generating naive summary statistics (NSS) for "+ (traitIdx[i] + 1) + "th variable in file '" + nssArgs.getPhenotypeFile() + "'.");
			nssResult = NewIt.newArrayList();
			naiveGWAS(pheFileIdx, traitIdx, i);
			printResult(traitIdx, i, pheFileIdx);
		}
		
		for (int i = 0; i < covIdx.length; i++)
		{
			Logger.printUserLog("");
			Logger.printUserLog("Generating naive summary statistics (NSS) for "+ (covIdx[i] + 1) + "th variable in file '" + nssArgs.getCovFile() + "'.");
			nssResult = NewIt.newArrayList();
			naiveGWAS(covFileIdx, covIdx, i);
			printResult(covIdx, i, covFileIdx);
		}

		String Fout = nssArgs.getOutRoot()  + ".list.nss";
		PrintStream nssList = FileUtil.CreatePrintStream(Fout);
		for (int i = 0; i < traitIdx.length; i++)
		{
			nssList.println(nssArgs.getOutRoot() + ".p." + (traitIdx[i] + 1) + ".nss");
		}
		for (int i = 0; i < covIdx.length; i++)
		{
			nssList.println(nssArgs.getOutRoot() + ".c." + (covIdx[i] + 1) + ".nss");
		}

		nssList.close();
		printCovMat();
	}

	private void printCovMat()
	{
		int[] pIdx = data.getMatchedSubjectIdx(pheFileIdx);
		int[] cIdx = data.getMatchedSubjectIdx(covFileIdx);

		double[][] pheVec = new double[data.getNumberOfSubjects()][traitIdx.length + covIdx.length];

		for (int subjectIdx = 0; subjectIdx < pIdx.length; subjectIdx++)
		{
			pheVec[subjectIdx][0]= this.data.getVariable(pheFileIdx, pIdx[subjectIdx], this.traitIdx[0]);
			for (int j = 0; j < covIdx.length; j++)
			{
				pheVec[subjectIdx][j+1] = this.data.getVariable(covFileIdx, cIdx[subjectIdx], this.covIdx[j]);
			}
		}

		PearsonsCorrelation pc; 
		RealMatrix mc;

		pc = new PearsonsCorrelation(pheVec);
		mc = pc.getCorrelationMatrix();

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

	private void naiveGWAS(int fileIdx, int[] variable,int tIdx) 
	{
		int monoLoci = 0;
		ChiSquaredDistributionImpl ci = new ChiSquaredDistributionImpl(1);
		ArrayList<SNP> snpList = this.mapFile.getMarkerList();

		int[] pIdx = data.getMatchedSubjectIdx(fileIdx);
		double[] Y = new double[pIdx.length];
		ArrayList<Double> pArray = NewIt.newArrayList();

		for (int subjectIdx = 0; subjectIdx < Y.length; subjectIdx++) 
		{
			Y[subjectIdx] = this.data.getVariable(fileIdx, pIdx[subjectIdx], variable[tIdx]);
		}
		Y = StatUtils.normalize(Y);

		int[] subIdx = data.getMatchedSubjectIdx(famFileIdx);
		for (int i = 0; i < this.gm.getNumMarker(); i++)
		{
			SNP snp = (SNP) snpList.get(i);
			
			if (frq[i][0] < nssArgs.getMAF() || frq[i][1] < nssArgs.getMAF() )
			{
				monoLoci++;
				continue;
			}

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
					int g = this.gm.getAdditiveScoreOnFirstAllele(subIdx[idx], i);
					if (g != 3) {
						sReg.addData(g, Y[j]);
						N += 1.0D;
						freq += g;
						mx += g;
						xx += g*g;
					}
				}
				double vg = (xx - (mx/N) * (mx/N) * N)/(N-1);

				if (vg == 0)
				{
					monoLoci++;
					continue;
				}
				if (N <= 1) 
				{
					monoLoci++;
					continue;
				}

				freq /= 2.0D * N;

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
			Logger.printUserLog("Removed " + monoLoci + " loci [MAF < " + nssArgs.getMAF() + "]." );			
		}
		else if (monoLoci == 1)
		{
			Logger.printUserLog("Removed " + monoLoci + " locus [MAF < " + nssArgs.getMAF() + "]." );			
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

	public void printResult(int[] variable, int tIdx, int fileIdx)
	{
		int fidx = variable[tIdx] + 1;
		String Fout=null;
		if (fileIdx == pheFileIdx)
		{
			Fout = nssArgs.getOutRoot()  + ".p." + fidx + ".nss";
		}
		else
		{
			Fout = nssArgs.getOutRoot()  + ".c." + fidx + ".nss";
		}
		PrintStream nssGWAS = FileUtil.CreatePrintStream(Fout);
		nssGWAS.println(OATHConst.SNP +"\t" + OATHConst.CHR + "\t" + OATHConst.BP+ "\t" + OATHConst.RefAle + "\t" + OATHConst.AltAle +"\t" + OATHConst.RAF + "\t" + OATHConst.Vg + "\t" + OATHConst.BETA +"\t" + OATHConst.SE + "\t" + OATHConst.CHI + "\t" + OATHConst.P);

		for (int i = 0; i < nssResult.size(); i++)
		{
			NSSGWASResult e1 = nssResult.get(i);
			nssGWAS.println(e1.printEGWASResult(lambdaGC));
		}
		nssGWAS.close();
		
		Logger.printUserLog("Write the NSS into '" + Fout + "'.");
	}
	
	public int getN()
	{
		return N;
	}
}
