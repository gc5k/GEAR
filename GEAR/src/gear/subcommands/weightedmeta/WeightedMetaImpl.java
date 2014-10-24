package gear.subcommands.weightedmeta;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;
import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;

import gear.gwassummary.GWASReader;
import gear.gwassummary.MetaStat;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.BufferedReader;
import gear.util.Logger;
import gear.util.NewIt;
import gear.util.SNPMatch;

public class WeightedMetaImpl extends CommandImpl
{

	@Override
	public void execute(CommandArguments cmdArgs)
	{
		wMetaArgs = (WeightedMetaArguments) cmdArgs;

		if (wMetaArgs.isQT())
		{
			Logger.printUserLog("Analysing summary statistics analysis for quantitative traits.\n");			
		}
		else
		{
			Logger.printUserLog("Analysing summary statistics analysis for case-contrl studies.\n");			
		}

		initial();

		if(wMetaArgs.getCMFile()!=null)
		{
			readCMFile();			
		}
		else
		{
			Logger.printUserLog("No correlation matrix is specified. The default correlation (unity matrix) will be used.");
		}
		MetaAnalysis();
	}

	private void initial()
	{
		gReader = new GWASReader(wMetaArgs.getMetaFile(), wMetaArgs.getKeys(), wMetaArgs.isQT(), wMetaArgs.isGZ());

		int NumMetaFile = gReader.getNumMetaFile();
		
		if (NumMetaFile < 2)
		{
			Logger.printUserError("At least two summary statistic files should be specified.\n");
			Logger.printUserError("GEAR quitted.\n");
			System.exit(0);
		}

		zMat = new double[NumMetaFile][NumMetaFile];
		for(int i = 0; i < zMat.length; i++)
		{
			zMat[i][i] = 1;
		}
	}

	private void readCMFile()
	{
		int NumMetaFile = gReader.getNumMetaFile();
		BufferedReader bf = BufferedReader.openTextFile(wMetaArgs.getCMFile(), "cm file.");
		Logger.printUserLog("Reading '" + wMetaArgs.getCMFile() + "'.");

		String[] d = null;
		int cnt=0;
		while ( (d = bf.readTokens())!= null ){
			if (d.length != NumMetaFile )
			{
				Logger.printUserError("incorrect '" + wMetaArgs.getCMFile() + "'.");
				System.exit(0);
			}
			for(int i = 0; i < d.length; i++)
			{
				zMat[cnt][i] = Double.parseDouble(d[i]);
			}
			cnt++;
		}	
	}

	private void MetaAnalysis()
	{
		int NumMetaFile = gReader.getNumMetaFile();

		double[][] Mx = new double[NumMetaFile][NumMetaFile];
		for (int i = 0; i < Mx.length; i++)
		{
			for (int j = 0; j < i + 1; j++)
			{
				Mx[i][j] = zMat[i][j];
				Mx[j][i] = zMat[i][j];
			}
		}

		RealMatrix gg = new Array2DRowRealMatrix(Mx);
		RealMatrix gg_Inv = (new LUDecompositionImpl(gg)).getSolver().getInverse();
		RealMatrix Unit = new Array2DRowRealMatrix(NumMetaFile, 1);
		for (int i = 0; i < Unit.getRowDimension(); i++)
		{
			Unit.setEntry(i, 0, 1);
		}
		RealMatrix tmp = Unit.transpose().multiply(gg_Inv);
		RealMatrix tmp1 = tmp.multiply(Unit);
		RealMatrix W = tmp.scalarMultiply(1/tmp1.getEntry(0, 0));

		Set<String> keys = gReader.getMetaSNPTable().keySet();
		for (Iterator<String> e=keys.iterator(); e.hasNext();)
		{
			String key = e.next();
			ArrayList<Integer> Int = gReader.getMetaSNPTable().get(key);
			if(Int.get(Int.size()-1) == NumMetaFile)
			{//common snp only
				GMRes gr = MetaCommon(key, Mx, W, Int);
				grArray.add(gr);
			}
			else
			{
				GMRes gr = MetaSNP(key, Mx, Int);
				grArray.add(gr);
			}
		}
		Collections.sort(grArray);
		PrintGMresults();
	}

	private GMRes MetaCommon(String key, double[][] Mx, RealMatrix W, ArrayList<Integer> Int)
	{
		StringBuffer direct = new StringBuffer();
		double gb = 0;
		double gse = 0;
		double[] se = new double[W.getColumnDimension()];
		GMRes gr = new GMRes(W.getColumnDimension());

		HashMap<String, MetaStat> m1 = gReader.getMetaStat().get(0);
		MetaStat ms1 = m1.get(key);
		gb += ms1.getEffect() * W.getEntry(0, 0);
		se[0] = ms1.getSE();
	
		gr.SetSNP(ms1.getSNP());
		gr.SetChr(ms1.getChr());
		gr.SetBP(ms1.getBP());
		gr.SetA1(ms1.getA1());
		gr.SetA2(ms1.getA2());
		if(ms1.getEffect() > 0)
		{
			direct.append("+");
		}
		else
		{
			direct.append("-");
		}
		for(int i = 1; i < gReader.getMetaStat().size(); i++)
		{
			HashMap<String, MetaStat> m2 = gReader.getMetaStat().get(i);
			MetaStat ms2 = m2.get(key); 
			se[i] = ms2.getSE();
			double b = 0;
			if (ms1.getA1() == ms2.getA1() || ms1.getA1() == SNPMatch.Flip(ms2.getA1())) //match A1 in the second meta
			{
				gb += ms2.getEffect() * W.getEntry(0, i);
				b = ms2.getEffect();
			}
			else if (ms1.getA1() == ms2.getA2() || ms1.getA1() == SNPMatch.Flip(ms2.getA2())) //match A2 in the second meta
			{
				gb += (-1) * ms2.getEffect() * W.getEntry(0, i);
				b = (-1) * ms2.getEffect();
			}
			if(ms2.getEffect() > 0)
			{
				direct.append("+");
			}
			else
			{
				direct.append("-");
			}
		}

		for(int i = 0; i < W.getRowDimension(); i++)
		{
			for(int j = 0; j < W.getRowDimension(); j++)
			{
				gse += W.getEntry(0, i) * W.getEntry(0, j) * se[i] * se[j] * Mx[i][j];
			}
		}
		gse = Math.sqrt(gse);
		double z= gb/gse;
		double p = 1;
		try
		{
			p=(1-unitNormal.cumulativeProbability(Math.abs(z)))*2;
		}
		catch (MathException e)
		{
			Logger.printUserError(e.toString());
		}
		gr.SetB(gb);
		gr.SetSE(gse);
		gr.SetP(p);
		gr.SetDirect(direct.toString());
		return gr;
	}

	private GMRes MetaSNP(String key, double[][] Mx, ArrayList<Integer> Int)
	{
		StringBuffer direct = new StringBuffer();
		for(int i = 0; i < Mx.length; i++)
		{
			direct.append("?");
		}
		GMRes gr = new GMRes(Int.get(Int.size()-1));
		int[] idx = new int[Int.get(Int.size()-1)];
		
		int cnt = 0;
		for(int i = 0; i < Int.size() - 1; i++)
		{
			if(Int.get(i) == 0) continue;
			idx[cnt++] = i;
		}

		double[][] mx = new double[idx.length][idx.length];
		for(int i = 0; i < idx.length; i++)
		{
			for(int j = 0; j < idx.length; j++)
			{
				mx[i][j] = Mx[idx[i]][idx[j]];
			}
		}

		RealMatrix gg = new Array2DRowRealMatrix(mx);
		RealMatrix gg_Inv = (new LUDecompositionImpl(gg)).getSolver().getInverse();
		RealMatrix Unit = new Array2DRowRealMatrix(mx.length, 1);
		for (int i = 0; i < Unit.getRowDimension(); i++)
		{
			Unit.setEntry(i, 0, 1);
		}
		RealMatrix tmp = Unit.transpose().multiply(gg_Inv);
		RealMatrix tmp1 = tmp.multiply(Unit);
		RealMatrix W = tmp.scalarMultiply(1/tmp1.getEntry(0, 0));

		double gb = 0;
		double gse = 0;
		double[] se = new double[W.getColumnDimension()];
		
		HashMap<String, MetaStat> m1 = gReader.getMetaStat().get(idx[0]);
		MetaStat ms1 = m1.get(key);
		gb += ms1.getEffect() * W.getEntry(0, 0);
		se[0] = ms1.getSE();

		gr.SetSNP(ms1.getSNP());
		gr.SetChr(ms1.getChr());
		gr.SetBP(ms1.getBP());
		gr.SetA1(ms1.getA1());
		gr.SetA2(ms1.getA2());
		if(ms1.getEffect() > 0)
		{
			direct.setCharAt(idx[0], '+');
		}
		else
		{
			direct.setCharAt(idx[0], '-');
		}

		for (int i = 1; i < idx.length; i++)
		{
			HashMap<String, MetaStat> m2 = gReader.getMetaStat().get(idx[i]);
			MetaStat ms2 = m2.get(key); 
			se[i] = ms2.getSE();
			double b = 0;
			if (ms1.getA1() == ms2.getA1() || ms1.getA1() == SNPMatch.Flip(ms2.getA1())) //match A1 in the second meta
			{
				gb += ms2.getEffect() * W.getEntry(0, i);
				b = ms2.getEffect();
			}
			else if (ms1.getA1() == ms2.getA2() || ms1.getA1() == SNPMatch.Flip(ms2.getA2())) //match A2 in the second meta
			{
				gb += (-1) * ms2.getEffect() * W.getEntry(0, i);
				b = -1 * ms2.getEffect();
			}
			if(ms1.getEffect() > 0)
			{
				direct.setCharAt(idx[i], '+');
			}
			else
			{
				direct.setCharAt(idx[i], '-');
			}

		}

		for (int i = 0; i < W.getRowDimension(); i++)
		{
			for (int j = 0; j < W.getRowDimension(); j++)
			{
				gse += W.getEntry(0, i) * W.getEntry(0, j) * se[i] * se[j] * Mx[i][j];
			}
		}
		gse = Math.sqrt(gse);
		double z= gb/gse;
		double p = 1;
		try
		{
			p=(1-unitNormal.cumulativeProbability(Math.abs(z)))*2;
		}
		catch (MathException e)
		{
			Logger.printUserError(e.toString());
		}
		gr.SetB(gb);
		gr.SetSE(gse);
		gr.SetZ(z);
		gr.SetP(p);
		gr.SetDirect(direct.toString());
		return gr;
	}

	private void PrintGMresults()
	{
        PrintWriter writer = null;
        try
        {
        	writer = new PrintWriter(new BufferedWriter(new FileWriter(wMetaArgs.getOutRoot()+".gmeta")));
        	Logger.printUserLog("Writting detailed test statistics into '"+wMetaArgs.getOutRoot() + ".gmeta.'\n");
        }
		catch (IOException e)
		{
			Logger.handleException(e, "An I/O exception occurred when writing '" + wMetaArgs.getOutRoot() + ".gmeta" + "'.\n");
		}

		for(int i = 0; i < grArray.size(); i++)
		{
			GMRes gr = grArray.get(i);
			if(i == 0)
			{
				writer.write(gr.printTitle() + "\n");
			}
			writer.write(gr.toString()+ "\n");
		}
		writer.close();
	}

	private WeightedMetaArguments wMetaArgs;
	private GWASReader gReader;
	private double[][] zMat;
	NormalDistributionImpl unitNormal = new NormalDistributionImpl(0, 1);
	private ArrayList<GMRes> grArray = NewIt.newArrayList();
}
