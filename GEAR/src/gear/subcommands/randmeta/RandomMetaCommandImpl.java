package gear.subcommands.randmeta;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ChiSquaredDistributionImpl;
import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math.stat.ranking.NaNStrategy;
import org.apache.commons.math.stat.ranking.NaturalRanking;
import org.apache.commons.math.stat.ranking.TiesStrategy;

import gear.gwassummary.GWASReader;
import gear.gwassummary.MetaStat;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.subcommands.lambdaD.LamUnit;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;
import gear.util.SNPMatch;

public class RandomMetaCommandImpl extends CommandImpl
{

	@Override
	public void execute(CommandArguments cmdArgs)
	{
		randArgs = (RandomMetaCommandArguments) cmdArgs;

		if (randArgs.isQT())
		{
			Logger.printUserLog("Analysing summary statistics analysis for quantitative traits.\n");			
		}
		else
		{
			Logger.printUserLog("Analysing summary statistics analysis for case-contrl studies.\n");			
		}

		initial();

//generating matrix
		String[] MetaFile = gReader.getMetaFile();
		for (int i = 0; i < MetaFile.length - 1; i++)
		{
			for (int j = (i + 1); j < MetaFile.length; j++)
			{
				Logger.printUserLog("File pair: " + (i + 1) + "-" + (j + 1));
				if (randArgs.isQT())
				{
					double[] size = randArgs.getQTsize();
					Kappa = 2 / (Math.sqrt(size[i] / size[j]) + Math
							.sqrt(size[j] / size[i]));
					Logger.printUserLog("Sample sizes for '" + MetaFile[i] + "': " + size[i]);
					Logger.printUserLog("Sample sizes for '" + MetaFile[j] + "': " + size[j]);
				}
				else
				{
					double[] size = randArgs.getCCsize();
					R1 = size[i * 2] / size[i * 2 + 1];
					R2 = size[j * 2] / size[j * 2 + 1];
					double s1 = size[i * 2] + size[i * 2 + 1];
					double s2 = size[j * 2] + size[j * 2 + 1];
					Kappa = 2 / (Math.sqrt(s1 / s2) + Math.sqrt(s2 / s1));
					Logger.printUserLog("Sample size for '" + MetaFile[i] + "': " + size[i * 2] + " cases, " + size[i * 2 + 1] + " controls; R1 = " + R1 + ".");
					Logger.printUserLog("Sample size for '" + MetaFile[j] + "': " + size[j * 2] + " cases, " + size[j * 2 + 1] + " controls; R2 = " + R2 + ".");
				}
				Logger.printUserLog("Kappa: " + Kappa);
				calculateLambdaD(i, j);
			}
		}
		WriteMat();
		Logger.printUserLog("=========================================================");
		Logger.printUserLog("Results has been saved in '" + randArgs
				.getOutRoot() + ".lmat'.");
		
	}
	
	private void initial()
	{
		boolean[] FileKeep = new boolean[randArgs.getMetaFile().length];
		Arrays.fill(FileKeep, true);
		gReader = new GWASReader(randArgs.getMetaFile(), FileKeep, randArgs.getKeys(), randArgs.isQT(), randArgs.isGZ(), randArgs.isChr(), randArgs.getChr());

		Me = randArgs.getMe();
		int NumMetaFile = randArgs.getMetaFile().length;

		if (NumMetaFile < 2)
		{
			Logger.printUserError("At least two summary statistic files should be specified.\n");
			Logger.printUserError("GEAR quitted.\n");
			System.exit(0);
		}

		bMat = new double[NumMetaFile][NumMetaFile];
		zscoreMat = new double[NumMetaFile][NumMetaFile];

		kMat = new double[NumMetaFile][NumMetaFile];

		for(int i = 0; i < NumMetaFile; i++)
		{
//			Arrays.fill(lamMat[i], 1);
			Arrays.fill(bMat[i], 1);
			Arrays.fill(zscoreMat[i], 1);
			Arrays.fill(kMat[i], 1);
		}

//reading meta files
	}

	private void calculateLambdaD(int idx1, int idx2)
	{
		ArrayList<LamUnit> LamArray = NewIt.newArrayList();
		BetaVec Bvec = new BetaVec();
//		ArrayList<BetaVec> BV = NewIt.newArrayList();
    	DescriptiveStatistics T0 = new DescriptiveStatistics();

		int cntAmbiguous = 0;
		HashMap<String, MetaStat> SumStat1 = gReader.getMetaStat().get(idx1);
		HashMap<String, MetaStat> SumStat2 = gReader.getMetaStat().get(idx2);

		ArrayList<String> snpArray = gReader.getMetaSNPArray().get(idx1);
		
		int[][] KeyIdx = gReader.getKeyIndex();
		for(String snp : snpArray)
		{
			if (!SumStat2.containsKey(snp) || !SumStat1.containsKey(snp))
			{
				continue;
			}
			MetaStat ms1 = SumStat1.get(snp);
			MetaStat ms2 = SumStat2.get(snp);
			double d = 0;
			
			if (KeyIdx[idx1][GWASReader.SE] != -1)
			{
				if (SNPMatch.isAmbiguous(ms1.getA1(), ms1.getA2()))
				{
					cntAmbiguous++;
					continue;
				}
			}
			if (KeyIdx[idx2][GWASReader.SE] != -1)
			{
				if (SNPMatch.isAmbiguous(ms2.getA1(), ms2.getA2()))
				{
					cntAmbiguous++;
					continue;
				}
			}

			if (ms1.getA1() == ms2.getA1() || ms1.getA1() == SNPMatch.Flip(ms2.getA1())) //match A1 in the second meta
			{
				d = (ms1.getEffect() - ms2.getEffect()) * (ms1.getEffect() - ms2.getEffect()) / (ms1.getSE() * ms1.getSE() + ms2.getSE() * ms2.getSE());
				Bvec.addStats(ms1.getEffect(), ms2.getEffect(), ms1.getSE(), ms2.getSE());
			}
			else if (ms1.getA1() == ms2.getA2() || ms1.getA1() == SNPMatch.Flip(ms2.getA2())) //match A2 in the second meta
			{
				d = ((-1) * ms1.getEffect() - ms2.getEffect()) * ((-1) * ms1.getEffect() - ms2.getEffect()) / (ms1.getSE() * ms1.getSE() + ms2.getSE() * ms2.getSE());
				Bvec.addStats(ms1.getEffect(), -1 * ms2.getEffect(), ms1.getSE(), ms2.getSE());
			}
			else
			{
				cntAmbiguous++;
				continue;
			}
			T0.addValue(d);
			LamArray.add(new LamUnit(d, ms1, ms2));
        }

		if (cntAmbiguous > 0)
		{
			if (cntAmbiguous == 1)
			{
				Logger.printUserLog("Removed " + cntAmbiguous + " ambiguous locus (AT/GC).");				
			}
			else
			{
				Logger.printUserLog("Removed " + cntAmbiguous + " ambiguous loci (AT/GC).");
			}
		}
		Logger.printUserLog("Found " + T0.getN() + " paired summary statistics between two files.");

//select independent snps
		double[] sortLD = T0.getSortedValues();
		double[] DesStat = null;
		int[] selIdx = null;

//		double[][] beta = null;
//		double[][] zscore = null;
		
		if (Me < 0)
		{//use all
			DesStat = new double[sortLD.length];
			System.arraycopy(sortLD, 0, DesStat, 0, sortLD.length);
			selIdx = new int[sortLD.length];
			for (int i = 0; i < sortLD.length; i++)
			{
				selIdx[i] = i;
			}
			Bvec.setSelected();

		}
		else if (sortLD.length <= Me)
		{//use available ones
			DesStat = new double[sortLD.length];
			System.arraycopy(sortLD, 0, DesStat, 0, sortLD.length);
			selIdx = new int[sortLD.length];
			for (int i = 0; i < sortLD.length; i++)
			{
				selIdx[i] = i;
			}
			Bvec.setSelected();

		}
		else
		{//use Me
			DesStat = new double[(int) Math.ceil(Me)];
			selIdx = new int[(int) Math.ceil(Me)];
			for (int i = 0; i < DesStat.length; i++)
			{
				selIdx[i] = (int) Math.floor( (i*1.0 + 1)/Me * sortLD.length ) -1;
				DesStat[i] = sortLD[selIdx[i]];
			}
			Bvec.setSelected(selIdx);
		}

		Bvec.CalCorrelation();
		Bvec.printOut();
		bMat[idx2][idx1] = Bvec.getBcorrelation();
		bMat[idx1][idx2] = Bvec.getPvBcorrelation();

		zscoreMat[idx2][idx1] = Bvec.getZcorrelation();
		zscoreMat[idx1][idx2] = Bvec.getPvZcorrelation();

		
		if (randArgs.isVerboseGZ())
		{
			VerboseGZ(LamArray, T0.getValues(), idx1, idx2);
		}
		else if (randArgs.isVerbose())
		{
			Verbose(LamArray, T0.getValues(), idx1, idx2);
		}
		else
		{
			NotVerbose(LamArray, T0.getValues(), idx1, idx2, selIdx);
		}
	}

	private void NotVerbose(ArrayList<LamUnit> LamArray, double[] ld, int idx1, int idx2, int[] selIdx)
	{
		Collections.sort(LamArray);
		ChiSquaredDistributionImpl chiDis = new ChiSquaredDistributionImpl(1);

        PrintWriter writer = null;
        try
        {
        	writer = new PrintWriter(new BufferedWriter(new FileWriter(randArgs.getOutRoot() + "." + (idx1+1) + "-" + (idx2+1) + ".rlam")));
        	Logger.printUserLog("Writting detailed test statistics into '"+randArgs.getOutRoot() + "." + (idx1+1) + "-" + (idx2+1) + ".rlam" + ".'\n");
        }
		catch (IOException e)
		{
			Logger.handleException(e, "An I/O exception occurred when writing '" + randArgs.getOutRoot() + "." + (idx1+1) + "-" + (idx2+1) + ".rlam" + "'.\n");
		}

        FileUtil.CreatePrintStream(new String(randArgs.getOutRoot() + "." + (idx1+1) + "-" + (idx2+1) + ".lam"));
       	writer.write(titleLine);

        for (int i = 0; i < selIdx.length; i++)
        {
        	LamUnit lu = LamArray.get(selIdx[i]);
        	MetaStat ms1 = lu.getMS1();
        	MetaStat ms2 = lu.getMS2();
        	double chi0 = 0;
			try
			{
				chi0 = chiDis.inverseCumulativeProbability((i+1)/(selIdx.length+0.05));;
			}
			catch (MathException e)
			{
				e.printStackTrace();
			}
			double lambda = lu.getChi1()/chi0;
        	writer.write(ms1.getSNP() + "\t" + ms1.getChr() + "\t" + ms1.getBP() + "\t" + ms1.getA1() + "\t" + ms1.getEffect() + "\t"  + ms1.getSE() + "\t" + ms1.getP() + "\t"+ ms2.getEffect() + "\t" + ms2.getSE() + "\t" + ms2.getP() + "\t" +lu.getChi1() + "\t" + chi0 + "\t" + lambda + "\n");
        }
        writer.close();
	}

	private void Verbose(ArrayList<LamUnit> LamArray, double[] ld, int idx1, int idx2)
	{
		ChiSquaredDistributionImpl chiDis = new ChiSquaredDistributionImpl(1);

		NaturalRanking ranking = new NaturalRanking(NaNStrategy.MINIMAL, TiesStrategy.SEQUENTIAL);
        double[] ranks = ranking.rank(ld);

        PrintWriter writer = null;
        try
        {
        	writer = new PrintWriter(new BufferedWriter(new FileWriter(randArgs.getOutRoot() + "." + (idx1+1) + "-" + (idx2+1) + ".rlam")));
        	Logger.printUserLog("Writting detailed test statistics into '" + randArgs.getOutRoot() + "." + (idx1+1) + "-" + (idx2+1) + ".rlam.'\n");
        }
		catch (IOException e)
		{
			Logger.handleException(e, "An I/O exception occurred when writing '" + randArgs.getOutRoot() + "." + (idx1+1) + "-" + (idx2+1) + ".rlam" + "'.");
		}

        FileUtil.CreatePrintStream(new String(randArgs.getOutRoot() + "." + (idx1+1) + "-" + (idx2+1) + ".lam"));
       	writer.write(titleLine);

        for (int i = 0; i < ranks.length; i++)
        {
        	LamUnit lu = LamArray.get(i);
        	MetaStat ms1 = lu.getMS1();
        	MetaStat ms2 = lu.getMS2();
        	double chi0 = 0;
			try
			{
				chi0 = chiDis.inverseCumulativeProbability(ranks[i]/(ranks.length+1));
			}
			catch (MathException e)
			{
				e.printStackTrace();
			}
			double lambda = lu.getChi1()/chi0;
        	writer.write(ms1.getSNP() + "\t" + ms1.getChr() + "\t" + ms1.getBP() + "\t" + ms1.getA1() + "\t" + ms1.getEffect() + "\t"  + ms1.getSE() + "\t" + ms1.getP() + "\t"+ ms2.getEffect() + "\t" + ms2.getSE() + "\t" + ms2.getP() + "\t" +lu.getChi1() + "\t" + chi0 + "\t" + lambda + "\n");
        }
        writer.close();
	}

	private void VerboseGZ(ArrayList<LamUnit> LamArray, double[] ld, int idx1, int idx2)
	{
		ChiSquaredDistributionImpl chiDis = new ChiSquaredDistributionImpl(1);

		NaturalRanking ranking = new NaturalRanking(NaNStrategy.MINIMAL, TiesStrategy.SEQUENTIAL);
        double[] ranks = ranking.rank(ld);
		BufferedWriter GZ = null;
		GZ = FileUtil.ZipFileWriter(new String(randArgs.getOutRoot() + "." + (idx1+1) + "-" + (idx2+1) + ".rlam.gz"));
    	Logger.printUserLog("Writting detailed test statistics into '"+randArgs.getOutRoot() + "." + (idx1+1) + "-" + (idx2+1) + ".rlam.gz.'\n");

       	try
		{
			GZ.append(titleLine);
		}
		catch (IOException e)
		{
			Logger.handleException(e, "error in writing " + new String(randArgs.getOutRoot() + "." + (idx1+1) + "-" + (idx2+1) + ".rlam.gz"));
		}

        for (int i = 0; i < ranks.length; i++)
        {
        	LamUnit lu = LamArray.get(i);
        	MetaStat ms1 = lu.getMS1();
        	MetaStat ms2 = lu.getMS2();
        	double chi0 = 0;
			try
			{
				chi0 = chiDis.inverseCumulativeProbability(ranks[i]/(ranks.length+1));
			}
			catch (MathException e)
			{
				e.printStackTrace();
			}
			double lambda = lu.getChi1()/chi0;
        	try
			{
				GZ.write(ms1.getSNP() + "\t" + ms1.getChr() + "\t" + ms1.getBP() + "\t" + ms1.getA1() + "\t" + ms1.getEffect() + "\t"  + ms1.getSE() + "\t" + ms1.getP() + "\t"+ ms2.getEffect() + "\t" + ms2.getSE() + "\t" + ms2.getP() + "\t" +lu.getChi1() + "\t" + chi0 + "\t" + lambda + "\n");
			}
			catch (IOException e)
			{
				Logger.handleException(e, "error in writing " + new String(randArgs.getOutRoot() + "." + (idx1+1) + "-" + (idx2+1) + ".lam.gz"));
			}
        }
        try
		{
			GZ.close();
		}
		catch (IOException e)
		{
			Logger.handleException(e, "error in writing " + new String(randArgs.getOutRoot() + "." + (idx1+1) + "-" + (idx2+1) + ".lam"));
		}
	}

	private void WriteMat()
	{
		//cm matrix
		PrintWriter cwriter = null;
		try
		{
			cwriter = new PrintWriter(new BufferedWriter(new FileWriter(randArgs.getOutRoot() + ".zcm")));
		}
		catch (IOException e)
		{
			Logger.handleException(e, "An I/O exception occurred when writing '" + randArgs.getOutRoot() + ".zcm" + "'.");
		}

		for (int i = 0; i < zscoreMat.length; i++)
		{
			for (int j = 0; j < zscoreMat[i].length; j++)
			{
				cwriter.print(String.format("%.4f", zscoreMat[i][j]) + " ");
			}
			cwriter.println();
		}
		cwriter.close();

		//Xmatrix
		PrintWriter bwriter = null;
		try
		{
			bwriter = new PrintWriter(new BufferedWriter(new FileWriter(randArgs.getOutRoot() + ".bcm")));
		}
		catch (IOException e)
		{
			Logger.handleException(e, "An I/O exception occurred when writing '" + randArgs.getOutRoot() + ".bcm" + "'.");
		}

		for (int i = 0; i < bMat.length; i++)
		{
			for (int j = 0; j < bMat[i].length; j++)
			{
				bwriter.print(String.format("%.4f", bMat[i][j]) + " ");
			}
			bwriter.println();
		}
		bwriter.close();
	}

	private double Me = 30000;

	private double R1 = 1;
	private double R2 = 1;
	private double Kappa = 1;
	private RandomMetaCommandArguments randArgs;

	private String titleLine= "SNP\tChr\tBp\tA1\tBETA1\tSE1\tP1\tBETA2\tSE2\tP2\tChiObs\tChiExp\tLambdaD\n";

	private GWASReader gReader;

	private double[][] bMat;
	private double[][] zscoreMat;
	private double[][] kMat;

}
