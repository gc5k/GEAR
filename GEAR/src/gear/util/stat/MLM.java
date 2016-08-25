package gear.util.stat;

import java.text.DecimalFormat;
import java.util.ArrayList;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.stat.StatUtils;

import gear.util.Logger;
import gear.util.NewIt;

public class MLM 
{
	public MLM(double[][] u2, double[] y, boolean isMINQUE)
	{
		A = new Array2DRowRealMatrix[2];
		A[0] = new Array2DRowRealMatrix(u2);
		A[1] = new Array2DRowRealMatrix(y.length, y.length);
		for(int i = 0; i < A[1].getColumnDimension(); i++)
		{
			A[1].setEntry(i, i, 1);
		}

		Y = new Array2DRowRealMatrix(y);

		double[][] x = new double[y.length][1];
		for(int i = 0; i < y.length; i++)
		{
			x[i][0] = 1;
		}
		X = new Array2DRowRealMatrix(x);

		this.isMINQUE = isMINQUE;
	}

	public MLM(double[][] u2, double[][] x, double[] y, boolean isMINQUE, int[] covIdx)
	{
		this.covIdx = new int[covIdx.length];
		System.arraycopy(covIdx, 0, this.covIdx, 0, covIdx.length);
		A = new Array2DRowRealMatrix[2];
		A[0] = new Array2DRowRealMatrix(u2);
		A[1] = new Array2DRowRealMatrix(y.length, y.length);
		for(int i = 0; i < A[1].getColumnDimension(); i++)
		{
			A[1].setEntry(i, i, 1);
		}

		Y = new Array2DRowRealMatrix(y);

		double[][] x_ = new double[x.length][x[0].length+1];
		for (int i = 0; i < x_.length; i++)
		{
			x_[i][0] = 1;
		}

		for (int i = 0; i < x.length; i++)
		{
			System.arraycopy(x[i], 0, x_[i], 1, x[i].length);
		}
		X = new Array2DRowRealMatrix(x_);

		this.isMINQUE = isMINQUE;
	}
	
	public MLM(double[][][] u3, double[] y, boolean isMINQUE)
	{
		A = new Array2DRowRealMatrix[u3.length+1];
		for(int i = 0; i < u3.length; i++)
		{
			A[i] = new Array2DRowRealMatrix(u3[i]);
		}
		A[A.length-1] = new Array2DRowRealMatrix(y.length, y.length);
		for(int i = 0; i < A[0].getColumnDimension(); i++)
		{
			A[A.length-1].setEntry(i, i, 1);
		}

		Y = new Array2DRowRealMatrix(y);

		double[][] x = new double[y.length][1];
		for(int i = 0; i < y.length; i++)
		{
			x[i][0] = 1;
		}
		X = new Array2DRowRealMatrix(x);
		this.isMINQUE = isMINQUE;

	}

	public MLM(double[][][] u3, double[][] x, double[] y, boolean isMINQUE, int[] covIdx)
	{
		this.covIdx = new int[covIdx.length];
		System.arraycopy(covIdx, 0, this.covIdx, 0, covIdx.length);

		A = new Array2DRowRealMatrix[u3.length+1];
		for (int i = 0; i < u3.length; i++)
		{
			A[i] = new Array2DRowRealMatrix(u3[i]);
		}
		A[A.length-1] = new Array2DRowRealMatrix(y.length, y.length);
		for (int i = 0; i < A[1].getColumnDimension(); i++)
		{
			A[A.length-1].setEntry(i, i, 1);
		}

		Y = new Array2DRowRealMatrix(y);

		double[][] x_ = new double[x.length][x[0].length+1];
		for (int i = 0; i < x_.length; i++)
		{
			x_[i][0] = 1;
		}

		for (int i = 0; i < x.length; i++)
		{
			System.arraycopy(x[i], 0, x_[i], 1, x[i].length);
		}
		X = new Array2DRowRealMatrix(x_);
		this.isMINQUE = isMINQUE;

	}

	public void MINQUE()
	{
		if (!isMINQUE)
		{
			Logger.printUserLog("REML estimation (may be slow if sample size is big...; negative VC is constrained to zero)");
		}
		else 
		{
			Logger.printUserLog("MINQUE estimation (may be slow if sample size is big...)");
		}

		double[] v = new double[A.length];
		double vt = StatUtils.variance(Y.getColumn(0));
		for(int i = 0; i < v.length; i++) v[i] = vt/v.length;

		Logger.printUserLog("Prior values for variance components are " + fmt.format(vt/v.length)+".");
		Logger.printUserLog("");

		String ts = new String("Iter.\tLog(L)\t");
		for(int i = 0; i < A.length-1; i++)
		{
			ts += "\t" + "V" + (i+1);
		}
		ts +="\tVe";
		
		if (!isMINQUE) Logger.printUserLog(ts);

		RealMatrix var = new Array2DRowRealMatrix(v);
		RealMatrix B = null;

		RealMatrix Int_V = null;
		RealMatrix X_IntV = null;
		RealMatrix X_IntV_X = null;
		RealMatrix X_IntV_X__Int = null;
		RealMatrix reml_m = null;

		int cnt = 0;
		while(true)
		{
			VAR.add(var.copy());
			
		//prepare V matrix
			RealMatrix V = null;
			for (int i = 0; i < var.getRowDimension(); i++)
			{
				if(i == 0) {
					V = A[i].scalarMultiply(var.getEntry(i, 0)).copy();
				}
				else {
					V=V.add(A[i].scalarMultiply(var.getEntry(i, 0)));
				}
			}
			Int_V =  (new LUDecompositionImpl(V)).getSolver().getInverse();
			X_IntV = (X.transpose()).multiply(Int_V);
			X_IntV_X = X_IntV.multiply(X);
			X_IntV_X__Int = (new LUDecompositionImpl(X_IntV_X)).getSolver().getInverse();

		//Estimate B using current value
			B = ((X_IntV_X__Int.multiply(X.transpose())).multiply(Int_V)).multiply(Y);
			RealMatrix yB = Y.subtract(X.multiply(B));
			RealMatrix re = (yB.transpose().multiply(Int_V)).multiply(yB);

		//Cal LOD, it's approximation, ignor determint(V)
			double lodProxi = -0.5 * (Y.getRowDimension() * Math.log(2*3.1415) + re.getEntry(0, 0));
			LOD.add(lodProxi);

			if (LOD.size() > 1)
			{
				double lod_diff = Math.abs(LOD.get(LOD.size() -1) - LOD.get(LOD.size()-2));
				if (lod_diff < converge) break;
			}

		//iteration
			RealMatrix tmp2 = (((Int_V.multiply(X)).multiply(X_IntV_X__Int)).multiply(X.transpose())).multiply(Int_V);
			RealMatrix Q = Int_V.subtract(tmp2);

			reml_m = new Array2DRowRealMatrix(A.length, A.length);
			RealMatrix reml_y = new Array2DRowRealMatrix(A.length, 1);

			for (int i = 0; i < reml_m.getRowDimension(); i++)
			{
				for (int j = 0; j < reml_m.getColumnDimension(); j++)
				{
					double ele = (((Q.multiply(A[j])).multiply(Q)).multiply(A[i])).getTrace();
					reml_m.setEntry(i, j, ele);
				}
				reml_y.setEntry(i, 0, (((((Y.transpose()).multiply(Q)).multiply(A[i])).multiply(Q)).multiply(Y)).getEntry(0, 0));
			}
			var = (reml_y.transpose().multiply((new LUDecompositionImpl(reml_m)).getSolver().getInverse())).transpose();
			
			for (int i = 0; i < var.getRowDimension(); i++)
			{
				if (var.getEntry(i, 0) < 0) var.setEntry(i, 0, 0);
			}

			cnt++;
			if (isMINQUE)
			{
				if (cnt == 2)
				{
					break;
				}
			}
			else
			{
				if (cnt == 1)
				{
					continue;
				}
				printIter(cnt);
			}
		}

		if (!isMINQUE) printIter(cnt+1);

		BETA = B;
		V_BETA = X_IntV_X__Int.multiply(X_IntV_X).multiply(X_IntV_X__Int);
		V_VAR = ((new LUDecompositionImpl(reml_m)).getSolver().getInverse()).scalarMultiply(2);
	}
	
	public void printIter(int cnt)
	{
		String its = new String();
		its +=  (new Integer(cnt-1)).toString() + "\t";
		its += fmt.format(LOD.get(LOD.size() -1));

		RealMatrix var_ = VAR.get(cnt-1);
		for (int i = 0; i < var_.getRowDimension(); i++)
		{
			its += "\t" + fmt.format(var_.getEntry(i, 0));
			if (var_.getEntry(i, 0) < 0)
			{
				var_.setEntry(i, 0, 0);
			}
		}
		Logger.printUserLog(its);

	}
	
	public void printVC()
	{
		//Summary VC
		Logger.printUserLog("");
		Logger.printUserLog("Summary result of VC analysis");
		Logger.printUserLog("----------------------------------------------------------------------------------------------");
		Logger.printUserLog("Comp.\tEstimate\tSE\tProp.");
		RealMatrix v_ = VAR.get(VAR.size()-1);

		double Vp = 0;
		for (int i = 0; i < v_.getRowDimension(); i++)
		{
			Vp += v_.getEntry(i, 0);
		}

		for (int i = 0; i < v_.getRowDimension() -1; i++)
		{
			Logger.printUserLog("V" + (i+1) + "\t" + fmt.format(v_.getEntry(i, 0)) + "\t" + fmt.format(Math.sqrt(V_VAR.getEntry(i, i)))+"\t" + fmt.format(v_.getEntry(i, 0)/Vp));
		}
		Logger.printUserLog("Ve" + "\t" + fmt.format(v_.getEntry(v_.getRowDimension()-1, 0))+ "\t" + fmt.format(Math.sqrt(V_VAR.getEntry(v_.getRowDimension()-1, v_.getRowDimension()-1))) +"\t" + fmt.format(v_.getEntry(v_.getRowDimension()-1, 0)/Vp) );
		Logger.printUserLog("Vp\t" + fmt.format(Vp));

		Logger.printUserLog("");
		Logger.printUserLog("Covariance structure for VC");

		for (int i = 0; i < V_VAR.getRowDimension(); i++)
		{
			String s = new String();
			for (int j = 0; j <=i; j++)
			{
				s += fmt.format(V_VAR.getEntry(i, j)) + "\t";
			}
			Logger.printUserLog(s);
		}
		Logger.printUserLog("----------------------------------------------------------------------------------------------");

		Logger.printUserLog("");
		Logger.printUserLog("Generalized linear square estimation (GLSE) for fixed effects");
		Logger.printUserLog("----------------------------------------------------------------------------------------------");
		Logger.printUserLog("Para.\tEstimate\tSE");

		for (int i = 0; i < BETA.getRowDimension(); i++)
		{
			if (i == 0)
			{
				Logger.printUserLog("Mean\t" + fmt.format(BETA.getEntry(i, 0)) + "\t" + fmt.format(Math.sqrt(V_BETA.getEntry(i, i))));
			}
			else 
			{
				Logger.printUserLog("Cov" + (covIdx[i-1]+1) +"\t"+ fmt.format(BETA.getEntry(i, 0)) + "\t" + fmt.format(Math.sqrt(V_BETA.getEntry(i, i))));
			}
		}
		Logger.printUserLog("");
		Logger.printUserLog("Covariance structrue for GLSE");
		for (int i = 0; i < V_BETA.getRowDimension(); i++)
		{
			String s = new String();
			for (int j = 0; j <= i; j++)
			{
				s += fmt.format(V_BETA.getEntry(i, j)) + "\t";
			}
			Logger.printUserLog(s);
		}
		Logger.printUserLog("----------------------------------------------------------------------------------------------");

	}

	public RealMatrix getVC()
	{
		return VAR.get(VAR.size() - 1);
	}

	public ArrayList<RealMatrix> getVCList()
	{
		return VAR;
	}

	public RealMatrix getVarVC()
	{
		return V_VAR;
	}

	public RealMatrix getBeta()
	{
		return BETA;
	}

	public RealMatrix getVarBeta()
	{
		return V_BETA;
	}

	private RealMatrix[] A = null;
	private RealMatrix Y = null;
	private RealMatrix X = null;
	private RealMatrix BETA = null;
	private RealMatrix V_BETA = null;
	private RealMatrix V_VAR = null;

	private ArrayList<Double> LOD = NewIt.newArrayList();
	private ArrayList<RealMatrix> VAR = NewIt.newArrayList();

	private int[] covIdx = null;
	private boolean isMINQUE = false;
	private double converge = 0.1;
	DecimalFormat fmt = new DecimalFormat("0.000");

}
