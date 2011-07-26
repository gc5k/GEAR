package score.dependonweka;
/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
import weka.core.matrix.Matrix;
/**
 * It implements AbstractScore.  It is developed to calculated score by logistic regression when without adjustment.
 * @author Guo-Bo Chen, Zhejiang University, gc5k@zju.edu.cn
 */
public class LogisticMeanScore implements AbstractScore {
	private double[] Residual;//residual score
	private Matrix Response;//Y

	/**
	 * Construct the object.
	 * @param Y Response.
	 */
	public LogisticMeanScore (double[][] Y)
	{
		Response = new Matrix(Y);
		Residual = new double[Response.getRowDimension()];
	}

	/**
	 * Get coefficients of the logistic regression model constructed.
	 */
	public double[] getCoefficients()
	{
		double mean[]= new double[Response.getColumnDimension()];
		double mu[] = new double[1];
		for( int i=0; i<Response.getColumnDimension(); i++)
		{
			for( int j=0; j<Response.getRowDimension(); j++)
			{
				mean[i]+=Response.get(j, i);
			}
			mean[i]/=Response.getRowDimension();
		}
		mu[0]=mean[0];
		return mu; 
	}

	/**
	 * Calculate score by a mean model.
	 */
	public void CalculateScore() {
		double mean[]= new double[Response.getColumnDimension()];
		for( int i=0; i<Response.getColumnDimension(); i++)
		{
			for( int j=0; j<Response.getRowDimension(); j++)
			{
				mean[i]+=Response.get(j, i);
			}
			mean[i]/=Response.getRowDimension();
			for( int j=0; j<Response.getRowDimension(); j++)
			{
				Residual[j] = Response.get(j, i)-mean[i];
			}
		}
	}

	/**
	 * Get score calculated by the mean model.
	 */
	public double[] GetScore() {
		return Residual;
	}

	/**
	 * Set ridge.
	 */
	public void SetRidge(double r) {
	}

	public static void main(String[] args)
	{
		double[][] response = { { 1 }, { 0 }, { 1 }, { 0 }, { 1 } };
		Matrix Y = new Matrix(response);
		LogisticMeanScore QS = new LogisticMeanScore(response);
		QS.CalculateScore();
		double[] res= QS.GetScore();
		for( int i=0; i<res.length; i++)
		{
			System.out.println(res[i]);
		}
	}
}
