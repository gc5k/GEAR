package gear.subcommands.glsmeta.util;

import gear.gwassummary.GWASReader;
import gear.gwassummary.MetaStat;
import gear.util.Logger;
import gear.util.SNPMatch;

import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;

public class MetaGLS
{

	public MetaGLS(String snp, ArrayList<Integer> Int, double[][] corMat,
			GWASReader gReader, boolean isGC, boolean isAdjOverlapping, double[][] covTable, boolean isCenter)
	{
		this.snp = snp;
		this.cohort = Int.get(Int.size() - 1).intValue();
		this.cohortIdx = new int[cohort];
		this.isGC = isGC;
		this.gc = new double[cohort];
		this.gReader = gReader;
		this.isCenter = isCenter;

		int cnt = 0;
		for (int i = 0; i < Int.size() - 1; i++)
		{
			if (Int.get(i) == 0) continue;
			cohortIdx[cnt++] = i;
		}

		setGC(this.isGC, gReader.GetGC());
		System.out.println(snp + ", cohort " + this.cohort);
		getY();
		System.out.print("Y: ");
		for(int i = 0; i < Y.length; i++) System.out.print(Y[i][0] + " ");
		System.out.println();

		double[][] covMat = setCovMat(corMat);
		for(int i = 0; i < covMat.length; i++)
		{
			for(int j = 0; j < covMat[i].length; j++)
			{
				System.out.print(covMat[i][j] + " ");
			}
			System.out.println();
		}

		RealMatrix gg = new Array2DRowRealMatrix(covMat);

		isNonSingular = (new LUDecompositionImpl(gg)).getSolver()
				.isNonSingular();

		if (isNonSingular)
		{
			Weight = getWeight(gg);
			System.out.print("W: ");
			for(int i = 0; i < Weight.length; i++)
			{
				System.out.print(Weight[i] + " ");
			}
			System.out.println();
		}

		if(covTable != null)
		{
			X = new double[cohortIdx.length][covTable[0].length + 1];
			for(int i = 0; i < X.length; i++)
			{
				X[i][0] = 1;
				for(int j = 0; j < covTable[i].length; j++)
				{
					X[i][1+j] = covTable[i][j];
				}
			}
		}
		else
		{
			if(!isGC)
			{
				X = new double[cohortIdx.length][1];
				for(int i = 0; i < X.length; i++) X[i][0] = 1;
			}
			else
			{
				X = new double[cohortIdx.length][2];
				for(int i = 0; i < X.length; i++)
				{
					X[i][0] = 1;
					X[i][1] = gc[cohortIdx[i]];
				}
			}
		}

		for(int i = 0; i < X.length; i++)
		{
			for(int j = 0; j < X[i].length; j++)
			{
				System.out.print(X[i][j] + " ");
			}
			System.out.println();
		}

		standardization();

		RealMatrix gg_1 = new Array2DRowRealMatrix(covMat);

		boolean isNonSingular_1 = (new LUDecompositionImpl(gg_1)).getSolver()
					.isNonSingular();

		if (isNonSingular_1)
		{
			RealMatrix gg_1_Inv = (new LUDecompositionImpl(gg_1)).getSolver().getInverse();

			RealMatrix y = new Array2DRowRealMatrix(Y);
			RealMatrix x = new Array2DRowRealMatrix(sX);
			RealMatrix X_T = x.transpose();

			RealMatrix X_T_RInv_X = (X_T.multiply(gg_1_Inv)).multiply(x);
			RealMatrix X_T_RInv_X__Inv = (new LUDecompositionImpl(X_T_RInv_X)).getSolver().getInverse();
			RealMatrix X_T_RInv_y = (X_T.multiply(gg_1_Inv)).multiply(y);

			B = X_T_RInv_X__Inv.multiply(X_T_RInv_y);
			Sigma = X_T_RInv_X__Inv;
			glsMetaRes.SetB(B);
			glsMetaRes.SetSE(Sigma);

		}

	}

	private void standardization()
	{
		double m_y = 0, sd_y = 0;
		double xp_y = 0;
		for(int i = 0; i < Y.length; i++)
		{
			m_y += Weight[i] * Y[i][0];
			xp_y += Weight[i] * Y[i][0] * Y[i][0];
		}
		sd_y = Math.sqrt(xp_y - m_y * m_y);
		System.out.println("Y: " + m_y + " " + xp_y + " " + sd_y);

//		adjFactor = new double[X[0].length];

		sX = new double[X.length][X[0].length];
		for(int i = 0; i < X.length; i++)
		{
			sX[i][0] = 1;
		}
		for(int i = 1; i < X[0].length; i++)
		{
			double m_x = 0, sd_x = 0;
			double xp_x = 0;

			for(int j = 0; j < X.length; j++)
			{
				m_x += Weight[j] * X[j][i];
				xp_x += Weight[j] * X[j][i] * X[j][i];
			}

			sd_x = Math.sqrt(xp_x - m_x * m_x);
			System.out.println("X: " + m_x + " " + xp_x + " " + sd_x);

			for(int j = 0; j < X.length; j++)
			{
				if (isCenter)
				{
					sX[j][i] = (X[j][i] - m_x) * sd_y/sd_x;					
				}
				else
				{
					sX[j][i] = (X[j][i] - m_x) * sd_y/sd_x + m_y;					
				}
			}
		}

		System.out.println("sX");

		for(int i = 0; i < sX.length; i++)
		{
			for(int j = 0; j < sX[i].length; j++)
			{
				System.out.print(sX[i][j] + " ");
			}
			System.out.println();
		}
	}

	private double[][] setCovMat(double[][] corMat)
	{
		covMat = new double[cohort][cohort];
		for (int i = 0; i < cohortIdx.length; i++)
		{
			MetaStat ms1 = gReader.getMetaStat().get(cohortIdx[i]).get(snp);
			for (int j = 0; j < cohortIdx.length; j++)
			{
				MetaStat ms2 = gReader.getMetaStat().get(cohortIdx[j]).get(snp);
				double adjcor = corMat[cohortIdx[i]][cohortIdx[j]];
				covMat[i][j] = adjcor * ms1.getSE() * ms2.getSE(); // * Math.sqrt(gc[i]) * Math.sqrt(gc[j]);
			}
		}
		return covMat;
	}

	private void getY()
	{
		Y = new double[cohortIdx.length][1];
		glsMetaRes = new GLSMetaRes(cohort);
		char sign = '+';

		direction = new StringBuffer();
		for(int i = 0; i < gReader.getCohortNum(); i++)
		{
			direction.append('?');
		}

		MetaStat ms = null;
		boolean isAmbiguousLocus = false;
		boolean match = true;
		for(int i = 0; i < cohortIdx.length; i++)
		{
			ms = gReader.getMetaStat().get(cohortIdx[i]).get(snp);
			if (i == 0)
			{
				Y[i][0] = ms.getEffect();
				glsMetaRes.SetSNP(ms.getSNP());
				glsMetaRes.SetChr(ms.getChr());
				glsMetaRes.SetBP(ms.getBP());
				glsMetaRes.SetA1(ms.getA1());
				glsMetaRes.SetA2(ms.getA2());
				isAmbiguousLocus = SNPMatch.isAmbiguous(ms.getA1(), ms.getA2());
			}
			else
			{
				match = SNPMatch.isAllelesMatchForTwoLoci(glsMetaRes.GetA1(), glsMetaRes.GetA2(), ms.getA1(), ms.getA2());
				if (!match)
				{
					match = SNPMatch.isAllelesFlipMatchForTwoLoci(glsMetaRes.GetA1(), glsMetaRes.GetA2(), ms.getA1(), ms.getA2());
				}
				
				if (match)
				{
					Y[i][0] = ms.getEffect();

					if (glsMetaRes.GetA1() == ms.getA1() || glsMetaRes.GetA1() == SNPMatch.Flip(ms.getA1())) //match A1 in the second meta
					{
					}
					else if (glsMetaRes.GetA1() == ms.getA2() || glsMetaRes.GetA1() == SNPMatch.Flip(ms.getA2())) //match A2 in the second meta
					{
						Y[i][0] = -1 * Y[i][0];
					}
				}
				else
				{
					sign = ',';
				}
			}

			if(match)
			{
				if(Y[i][0] == 0)
				{
					sign = '0';
				}
				else if(Y[i][0] > 0)
				{
					sign = '+';
				}
				else
				{
					sign = '-';
				}
			}
			direction.setCharAt(cohortIdx[i], sign);
		}
		System.out.println(direction);
		glsMetaRes.SetAmbi(isAmbiguousLocus);
		glsMetaRes.SetDirect(direction.toString());

	}
	
	private double[] getWeight(RealMatrix gg)
	{
		RealMatrix gg_Inv = (new LUDecompositionImpl(gg)).getSolver()
				.getInverse();
		RealMatrix Unit = new Array2DRowRealMatrix(covMat.length, 1);
		for (int i = 0; i < Unit.getRowDimension(); i++)
		{
			Unit.setEntry(i, 0, 1);
		}
		RealMatrix tmp = Unit.transpose().multiply(gg_Inv);
		RealMatrix tmp1 = tmp.multiply(Unit);

		RealMatrix W = tmp.scalarMultiply(1 / tmp1.getEntry(0, 0));

		double gse = 1 / tmp1.getEntry(0, 0);
		if (gse < 0)
		{
			Logger.printUserLog("This locus has negative variance: " + gse);
		}
		gse = Math.sqrt(gse);
		return W.getRow(0);
	}

	private void setGC(boolean isGC, double[] Egc)
	{
		if (!isGC)
		{
			Arrays.fill(this.gc, 1);
		}
		else
		{
			for (int i = 0; i < cohortIdx.length; i++)
			{
				gc[i] = Egc[cohortIdx[i]];
			}
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

	public GLSMetaRes getGLSMetaRes()
	{
		return glsMetaRes;
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
	private GWASReader gReader;
	private boolean isCenter;
	
	private double[][] Y;
	private double[][] X;
//	private double[] adjFactor;
	private double[][] sX;
	
	private RealMatrix B;
	private RealMatrix Sigma;
	StringBuffer direction;
	
	private GLSMetaRes glsMetaRes;

}
