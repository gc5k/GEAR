package gear.subcommands.weightedmeta;

import java.util.ArrayList;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;

import gear.gwassummary.GWASReader;
import gear.gwassummary.MetaStat;

public class CovMatrix
{
	public CovMatrix(String snp, ArrayList<Integer> Int, double[][] corMat, GWASReader gReader)
	{
		this.snp = snp;
		this.cohort = Int.get(Int.size() -1).intValue();
		this.cohortIdx = new int[cohort];

		int cnt = 0;
		for(int i = 0; i < Int.size()-1; i++)
		{
			if(Int.get(i) == 0) continue;
			cohortIdx[cnt++] = i;
		}

		double[][] covMat = new double[cohort][cohort];
		for(int i = 0; i < cohortIdx.length; i++)
		{
			MetaStat ms1 = gReader.getMetaStat().get(cohortIdx[i]).get(snp);
			for(int j = 0; j < cohortIdx.length; j++)
			{
				MetaStat ms2 = gReader.getMetaStat().get(cohortIdx[j]).get(snp);
				covMat[i][j] = corMat[cohortIdx[i]][cohortIdx[j]] * ms1.getSE() * ms2.getSE();
			}
		}

		RealMatrix gg = new Array2DRowRealMatrix(covMat);
		RealMatrix gg_Inv = (new LUDecompositionImpl(gg)).getSolver().getInverse();
		RealMatrix Unit = new Array2DRowRealMatrix(covMat.length, 1);
		for(int i = 0; i < Unit.getRowDimension(); i++)
		{
			Unit.setEntry(i, 0, 1);
		}
		RealMatrix tmp = Unit.transpose().multiply(gg_Inv);
		RealMatrix tmp1 = tmp.multiply(Unit);
		RealMatrix W = tmp.scalarMultiply(1/tmp1.getEntry(0, 0));

		for (int i = 0; i < W.getRowDimension(); i++)
		{
			for (int j = 0; j < W.getRowDimension(); j++)
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

	private String snp;
	private double[][] covMat;
	private double gse;
	private int cohort;
	private double[] Weight;
	private int[] cohortIdx;
}
