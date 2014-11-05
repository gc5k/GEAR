package gear.subcommands.weightedmeta;

import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.EigenDecompositionImpl;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;

import gear.gwassummary.GWASReader;
import gear.gwassummary.MetaStat;

public class CovMatrix
{
	public CovMatrix(String snp, ArrayList<Integer> Int, double[][] corMat, GWASReader gReader, boolean isGC, boolean isGCInflationOnly)
	{
		this.snp = snp;
		this.cohort = Int.get(Int.size() -1).intValue();
		this.cohortIdx = new int[cohort];
		this.isGC = isGC;
		this.isGCInflationOnly = isGCInflationOnly;
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
				gc[i] = Egc[cohortIdx[i]];
				if(this.isGCInflationOnly)
				{
					gc[i] = gc[i] > 1 ? gc[i]:1;
				}				
			}
		}
		else
		{
			Arrays.fill(this.gc, 1);
		}

		double[][] covMat = new double[cohort][cohort];
		for(int i = 0; i < cohortIdx.length; i++)
		{
			MetaStat ms1 = gReader.getMetaStat().get(cohortIdx[i]).get(snp);
			for(int j = 0; j < cohortIdx.length; j++)
			{
				MetaStat ms2 = gReader.getMetaStat().get(cohortIdx[j]).get(snp);
				covMat[i][j] = corMat[cohortIdx[i]][cohortIdx[j]] * ms1.getSE() * ms2.getSE() * Math.sqrt(gc[i]) *  Math.sqrt(gc[j]);
			}
		}

		RealMatrix gg = new Array2DRowRealMatrix(covMat);

		double[] eigen = (new EigenDecompositionImpl(gg, 0.00000001)).getRealEigenvalues();
		for(int i = 0; i < eigen.length; i++)
		{
			System.out.println(eigen[i]);
		}

		RealMatrix gg_Inv = (new LUDecompositionImpl(gg)).getSolver().getInverse();
		RealMatrix Unit = new Array2DRowRealMatrix(covMat.length, 1);
		for(int i = 0; i < Unit.getRowDimension(); i++)
		{
			Unit.setEntry(i, 0, 1);
		}
		RealMatrix tmp = Unit.transpose().multiply(gg_Inv);
		RealMatrix tmp1 = tmp.multiply(Unit);
		RealMatrix W = tmp.scalarMultiply(1/tmp1.getEntry(0, 0));

		for (int i = 0; i < covMat.length; i++)
		{
			for (int j = 0; j < covMat[i].length; j++)
			{
				gse += W.getEntry(0, i) * W.getEntry(0, j) * covMat[i][j];
			}
		}
		gse = Math.sqrt(gse);
		Weight = W.getRow(0);
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

	private boolean isGC;
	private boolean isGCInflationOnly;
	private double[] gc;
	private String snp;
	private double[][] covMat;
	private double gse;
	private int cohort;
	private double[] Weight;
	private int[] cohortIdx;
}
