package gear.subcommands.mlmmeta.util;

import java.util.ArrayList;
import java.util.Arrays;

import gear.gwassummary.GWASReader;
import gear.gwassummary.MetaStat;

public class MLMCovMatrix
{
	public MLMCovMatrix(String snp, ArrayList<Integer> Int, double[] qtSize, double[][]mlmMat, double[][] corMat, GWASReader gReader, boolean isGC, boolean isGCALL, boolean isAdjOverlapping)
	{
		this.snp = snp;
		this.cohort = Int.get(Int.size() -1).intValue();
		this.cohortIdx = new int[cohort];
		this.isGC = isGC;
		this.isGCALL = isGCALL;
		this.gc = new double[cohort];
		
		int cnt = 0;
		for(int i = 0; i < Int.size()-1; i++)
		{
			if(Int.get(i) == 0) continue;
			cohortIdx[cnt++] = i;
		}

		if (this.isGC)
		{
			double[] Egc = gReader.GetGC();
			for(int i = 0; i < cohortIdx.length; i++)
			{
				if (this.isGCALL)
				{
					gc[i] = Egc[cohortIdx[i]];
				}
				else
				{
					gc[i] = Egc[cohortIdx[i]] > 1 ? Egc[cohortIdx[i]]:1;					
				}
			}
		}
		else
		{
			Arrays.fill(this.gc, 1);
		}

		covMat = new double[cohort][cohort];
		mMat = new double[cohort][cohort];
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
				mMat[i][j] = mlmMat[cohortIdx[i]][cohortIdx[j]] * Math.sqrt(1/qtSize[cohortIdx[i]]) * Math.sqrt(1/qtSize[cohortIdx[j]]);
			}
		}
	}

	public double[][] getA()
	{
		return mMat;
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
	private boolean isGCALL;
	private double[] gc;
	private double[] qtSize;
	private String snp;
	private double[][] mMat;
	private double[][] covMat;
	private double gse = 0;
	private int cohort;
	private double[] Weight;
	private int[] cohortIdx;
	private boolean isNonSingular;
}
