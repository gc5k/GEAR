package gear.subcommands.oath.synthesize.freader;

import java.util.ArrayList;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;

public class SynthMatrix
{
	public SynthMatrix(int N, String snp, ArrayList<Integer> Int, double[][] corMat, SynthFReader fReader)
	{
		this.N = N;
		this.snp = snp;
		this.cohort = Int.get(Int.size() -1).intValue();
		this.cohortIdx = new int[cohort];
		int cnt = 0;
		for(int i = 0; i < Int.size()-1; i++)
		{
			if(Int.get(i) == 0) continue;
			cohortIdx[cnt++] = i;
		}

		initial(Int, corMat, fReader);
		cal();
	}

	private void cal()
	{
		RealMatrix lMat = new Array2DRowRealMatrix(LambdaMat);
		RealMatrix oMat = new Array2DRowRealMatrix(OmegaMat);
		RealMatrix bMat = new Array2DRowRealMatrix(BETA);

		isNonSingular = (new LUDecompositionImpl(lMat)).getSolver().isNonSingular();
		RealMatrix lMat_Inv = (new LUDecompositionImpl(lMat)).getSolver().getInverse();

		RealMatrix tmp = oMat.multiply(bMat);

		OATHB = lMat_Inv.multiply(tmp);

		RealMatrix OATHB_t = OATHB.transpose();
		OATHBse = Math.sqrt(lMat_Inv.scalarMultiply( (1 - (OATHB_t.multiply(tmp)).getEntry(0, 0))/(N - lMat.getColumnDimension() ) ).getEntry(0, 0));
	}

	private void initial(ArrayList<Integer> Int, double[][] corMat, SynthFReader fReader)
	{
		LambdaMat = new double[cohort][cohort];

		for (int i = 0; i < cohortIdx.length; i++)
		{
			for (int j = 0; j < cohortIdx.length; j++)
			{
				LambdaMat[i][j] = corMat[cohortIdx[i]][cohortIdx[j]];
			}
		}

		OmegaMat = new double[cohort][cohort];
		BETA = new double[cohort][1];
		for (int i = 0; i < cohortIdx.length; i++)
		{
			SynthFStat ms1 = fReader.getMetaStat().get(cohortIdx[i]).get(snp);
			if (i == 0)
			{
				BETA[i][0] = ms1.getEffect();
			}
			else
			{
				BETA[i][0] = corMat[cohortIdx[i]][0];
			}
		}

		for (int i = 0; i < cohortIdx.length; i++)
		{
			SynthFStat ms1 = fReader.getMetaStat().get(cohortIdx[i]).get(snp);
			if ( i == 0)
			{
				LambdaMat[i][i] = ms1.getVG();
			}
			else
			{
				LambdaMat[0][i] = LambdaMat[i][0] = ms1.getEffect() * LambdaMat[0][0];
			}			
			OmegaMat[i][i] = LambdaMat[i][i];
		}
	}

	public double[][] getLambdaMatrix()
	{
		return LambdaMat;
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

	public double getOATHB()
	{
		return OATHB.getEntry(0, 0);
	}
	
	public double getOATHse()
	{
		return OATHBse;
	}

	private String snp;
	private double[][] LambdaMat;
	private double[][] OmegaMat;
	private double[][] BETA;
	private double gse = 0;
	private int cohort;
	private double[] Weight;
	private int[] cohortIdx;
	private int N;
	private boolean isNonSingular;
	private RealMatrix OATHB;
	private double OATHBse;
}
