package gear.util.stat;

import java.text.DecimalFormat;
import java.util.ArrayList;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.EigenDecomposition;
import org.apache.commons.math.linear.EigenDecompositionImpl;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;

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
	}

	public MLM(double[][][] u3, double[][] x, double[] y, boolean isMINQUE)
	{
		A = new Array2DRowRealMatrix[u3.length+1];
		for (int i = 0; i < u3.length; i++)
		{
			A[i] = new Array2DRowRealMatrix(u3[i]);
		}
		A[A.length] = new Array2DRowRealMatrix(y.length, y.length);
		for (int i = 0; i < A[1].getColumnDimension(); i++)
		{
			A[A.length].setEntry(i, i, 1);
		}

		Y = new Array2DRowRealMatrix(this.y);

		double[][] x_ = new double[x.length+1][x[0].length];
		for (int i = 0; i < x_.length; i++)
		{
			x_[i][0] = 1;
		}

		for (int i = 0; i < x.length; i++)
		{
			System.arraycopy(x[i], 0, x_[i+1], 1, x[i].length);
		}
		X = new Array2DRowRealMatrix(x_);
	}

	public void MINQUE()
	{
		double[] v = new double[A.length];
		for(int i = 0; i < v.length; i++) v[i] = 1.0/v.length;

		RealMatrix var = new Array2DRowRealMatrix(v);
		RealMatrix B = null;

		ArrayList<Double> LOD = NewIt.newArrayList();
		double lod_diff = 1;

		int cnt = 0;
		while(true)
		{
			if (! isMINQUE)
			{
				Logger.printUserLog("REML Iteration " + (cnt+1));
			}
			else 
			{
				Logger.printUserLog("MINQUE procedure");
			}

			VAR.add(var.copy());
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
//			System.out.println(V);

			RealMatrix Int_V =  (new LUDecompositionImpl(V)).getSolver().getInverse();
			RealMatrix t1 = (X.transpose()).multiply(Int_V);
			RealMatrix t2 = t1.multiply(X);
			RealMatrix tmp = (new LUDecompositionImpl(t2)).getSolver().getInverse();

		//Cal LOD
		//Estimate B
			B = ((tmp.multiply(X.transpose())).multiply(Int_V)).multiply(Y);
			RealMatrix yB = Y.subtract(X.multiply(B));
			RealMatrix re = (yB.transpose().multiply(Int_V)).multiply(yB);
//			double v_det = (new LUDecompositionImpl(V)).getDeterminant();
//			double v_det = getSumEigenValues(V);
			double lodProxi = -0.5 * (Y.getRowDimension() * Math.log(2*3.1415) + re.getEntry(0, 0));
			LOD.add(lodProxi);
			if (isMINQUE && cnt == 1) break;
			if (LOD.size() > 1)
			{
				lod_diff = Math.abs(LOD.get(LOD.size() -1) - LOD.get(LOD.size()-2));
				if(lod_diff < 0.1) break;
			}

		//MINQUE
			RealMatrix tmp2 = (((Int_V.multiply(X)).multiply(tmp)).multiply(X.transpose())).multiply(Int_V);
			RealMatrix Q = Int_V.subtract(tmp2);

			RealMatrix reml_m = new Array2DRowRealMatrix(A.length, A.length);
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

			///
			Logger.printUserLog(var.toString());

			for (int i = 0; i < var.getRowDimension(); i++)
			{
				if (var.getEntry(i, 0) < 0)
				{
					Logger.printUserLog(i+1 + " VC was constrained to 0.");
					var.setEntry(i, 0, 0);
				}
			}

			cnt++;
		}

		String s = new String("\n");
		s += "Iteration\tLog(L)";
		for (int i = 0; i <var.getRowDimension()-1; i++)
		{
			s += "\t" + "V" + (i+1);
		}
		s += "\tVe\n";
		DecimalFormat fmt = new DecimalFormat("0.000");

		for (int i = 1; i < LOD.size(); i++)
		{
			s +=i +"\t" + fmt.format(LOD.get(i)) + "\t";
			RealMatrix v_ = VAR.get(i);
			for(int j = 0; j < v_.getRowDimension(); j++)
			{
				s += fmt.format(v_.getEntry(j, 0)) + " ";
			}
			s += "\n";
		}
		Logger.printUserLog(s);
	}

	public RealMatrix getVar()
	{
		return VAR.get(VAR.size() - 1);
	}
	
	private double getSumEigenValues(RealMatrix V)
	{//det = sum(EigenValue)
		EigenDecomposition Ed = new EigenDecompositionImpl(V, 1e-5);
		double[] ev = Ed.getRealEigenvalues();
		double v_det = 0;

		for(int i = 0; i < ev.length; i++)
		{
			if(ev[i] > 0)
			{
				v_det += ev[i];
			}
			else
			{
				break;
			}
		}
		return v_det;
	}
	
	private double[] y = null;
	private RealMatrix[] A = null;
	private RealMatrix Y = null;
	private RealMatrix X = null;
	ArrayList<RealMatrix> VAR = NewIt.newArrayList();

	private boolean isMINQUE = false;

}
