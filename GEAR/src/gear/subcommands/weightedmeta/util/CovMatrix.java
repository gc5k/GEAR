package gear.subcommands.weightedmeta.util;

import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;

import gear.gwassummary.GWASReader;
import gear.gwassummary.MetaStat;
import gear.util.Logger;

public class CovMatrix
{
	public CovMatrix(String snp, ArrayList<Integer> Int, double[][] corMat, GWASReader gReader, boolean isGC, boolean isAdjOverlapping)
	{
		this.snp = snp;
		this.cohort = Int.get(Int.size() -1).intValue();
		this.cohortIdx = new int[cohort];
		this.isGC = isGC;
		this.gc = new double[cohort];
		
		int cnt = 0;
		for(int i = 0; i < Int.size()-1; i++)
		{
			if(Int.get(i) == 0) continue;
			cohortIdx[cnt++] = i;
		}

		if(this.isGC)
		{
			double[] Egc = gReader.GetGC();
			for(int i = 0; i < cohortIdx.length; i++)
			{
				gc[i] = Egc[cohortIdx[i]] > 1 ? Egc[cohortIdx[i]]:1;
			}
		}
		else
		{
			Arrays.fill(this.gc, 1);
		}

		covMat = new double[cohort][cohort];
		for(int i = 0; i < cohortIdx.length; i++)
		{
			MetaStat ms1 = gReader.getMetaStat().get(cohortIdx[i]).get(snp);
			for(int j = 0; j < cohortIdx.length; j++)
			{
				MetaStat ms2 = gReader.getMetaStat().get(cohortIdx[j]).get(snp);
				double adjcor = corMat[cohortIdx[i]][cohortIdx[j]];
				if (isAdjOverlapping)
				{
					adjcor = corMat[cohortIdx[i]][cohortIdx[j]] < 0 ? 0:adjcor;
				}
				covMat[i][j] = adjcor * ms1.getSE() * ms2.getSE() * Math.sqrt(gc[i]) *  Math.sqrt(gc[j]);
			}
		}
		
		RealMatrix gg = new Array2DRowRealMatrix(covMat);

		isNonSingular = (new LUDecompositionImpl(gg)).getSolver().isNonSingular();
		
		if(isNonSingular)
		{
			RealMatrix gg_Inv = (new LUDecompositionImpl(gg)).getSolver().getInverse();
			RealMatrix Unit = new Array2DRowRealMatrix(covMat.length, 1);
			for(int i = 0; i < Unit.getRowDimension(); i++)
			{
				Unit.setEntry(i, 0, 1);
			}
			RealMatrix tmp = Unit.transpose().multiply(gg_Inv);
			RealMatrix tmp1 = tmp.multiply(Unit);
			
			RealMatrix W = tmp.scalarMultiply(1/tmp1.getEntry(0, 0));
			
			double gse = 1/tmp1.getEntry(0, 0);
			if(gse < 0)
			{
				Logger.printUserLog("This locus has negative variance: " + gse);
			}
			gse = Math.sqrt(gse);
			Weight = W.getRow(0);			
		}
	}

	public double[][] getCovMatrix()
	{
		return covMat;
	}

	public double[] getWeights()
	{
		return Weight;
	}

	public double getGSE()
	{
		return gse;
	}

	public int[] getCohortIdx()
	{
		return cohortIdx;
	}

	public String getSNP()
	{
		return snp;
	}

	public boolean isNonSingular()
	{
		return isNonSingular;
	}

	private boolean isGC;
	private double[] gc;
	private String snp;
	private double[][] covMat;
	private double gse = 0;
	private int cohort;
	private double[] Weight;
	private int[] cohortIdx;
	private boolean isNonSingular;
}
