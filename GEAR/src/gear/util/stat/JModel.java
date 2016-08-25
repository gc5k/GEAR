package gear.util.stat;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;

import gear.util.Logger;

public class JModel 
{
	public JModel(double[][] XX, double[] b)
	{
		if (b.length > 1)
		{
			RealMatrix Mat_XX = new Array2DRowRealMatrix(XX);
			if ( (new LUDecompositionImpl (Mat_XX).getSolver().isNonSingular()) )
			{
				RealMatrix Dxx = new Array2DRowRealMatrix(Mat_XX.getRowDimension(), Mat_XX.getColumnDimension());
				for(int i = 0; i < Dxx.getColumnDimension(); i++)
				{
					Dxx.setEntry(i, i, Mat_XX.getEntry(i, i));
				}
				RealMatrix Mat_XX_Inv = (new LUDecompositionImpl(Mat_XX)).getSolver()
						.getInverse();

				RealMatrix Mat_B = new Array2DRowRealMatrix(b);
				JB = Mat_XX_Inv.multiply(Dxx).multiply(Mat_B);				
			}
			else
			{
				isSingular = true;
				Logger.printUserLog("Matrix is singular. Joint model is cancelled.");
				JB = new Array2DRowRealMatrix(b);				
			}

		}
		else
		{
			JB = new Array2DRowRealMatrix(b);
		}
	}

	public boolean isSingular()
	{
		return isSingular;
	}

	public double[] getJB()
	{
		return JB.getColumn(0);
	}

	private RealMatrix JB = null;
	private boolean isSingular = false;
}
