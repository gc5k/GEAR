package score;

import java.lang.Integer;
import java.lang.Math;
import weka.core.*;
import weka.core.matrix.Matrix;
import weka.classifiers.functions.Logistic;

/**
 * It implement AbstractScore. It is utilized when logistic regression should be used to build score. 
 * @author Guo-Bo Chen, Zhejiang University, gc5k@zju.edu.cn
 */
public class LogisticScore implements AbstractScore{
	private double[] Residual;// residual score
	private Matrix Response;// Y
	private Matrix Observe;// X
	private Attribute[] AttriX;// Attribute X
	private Attribute AttriY;// Attribute Y
	private Instances LS;// dataset
	private double ridge;// ridge, default 0;
	boolean intercept;// intercept, default true;

	/**
	 * Construct the object, and initialize the predictors and response.
	 * @param X Predictors
	 * @param Y Response
	 */
	public LogisticScore(double[][] X, double[][] Y) {
		ridge = 0;
		Observe = new Matrix(X);
		Response = new Matrix(Y);
		Residual = new double[Response.getRowDimension()];
		ridge = 0;
		intercept = true;
		CreateInstance();
	}

	/**
	 * Create instance with the given predictors and response.
	 */
	private void CreateInstance() {
		AttriX = new Attribute[Observe.getColumnDimension()];
		FastVector my_nominal_values = new FastVector(2);
		my_nominal_values.addElement("1");// 1 for case
		my_nominal_values.addElement("0");// 0 for control
		Attribute AttriY = new Attribute("Response", my_nominal_values);
		FastVector attributes = new FastVector(Observe.getColumnDimension()
				+ Response.getColumnDimension());
		for (int i = 0; i < AttriX.length; i++) {
			AttriX[i] = new Attribute("X" + (String) (new Integer(i).toString()));
			attributes.addElement(AttriX[i]);
		}
		attributes.addElement(AttriY);
		LS = new Instances("LS", attributes, Observe.getRowDimension());
		LS.setClassIndex(AttriY.index());

		for (int i = 0; i < Observe.getRowDimension(); i++) {
			Instance ins = new Instance(Observe.getColumnDimension()
					+ Response.getColumnDimension());
			for (int j = 0; j < Observe.getColumnDimension(); j++) {
				ins.setValue(AttriX[j], Observe.get(i, j));
			}
			ins.setValue(AttriY, Response.get(i, 0));
			ins.setDataset(LS);
			LS.add(ins);
		}
	}

	/**
	 * Set ridge of the logistic regression model.
	 */
	public void SetRidge(double r) {
		ridge = r;
	}

	/**
	 * Get score calculated by the logistic regression model.
	 */
	public double[] GetScore() {
		return Residual;
	}

	/**
	 * Get the coefficients of the logistic regression model constructed.
	 */
	public double[] getCoefficients()
	{
		Logistic Reg = new Logistic();
		try {
			Reg.buildClassifier(LS);
		} catch (Exception err) {
			err.printStackTrace();
		}
		// add by magic
		String s = Reg.toString();
		String[] ss = s.split("\n");
		
		String[] coeffi = new String[LS.numAttributes()+1-1];
		for (int i = 0; i < coeffi.length; i++) {
			coeffi[i] = ss[i + 3].trim();
			coeffi[i] = coeffi[i].split("\\s+")[1];
		}
		double[] coe = stringTOdouble(coeffi);
		for( int i=0; i<coe.length; i++)
		{
			coe[i]*=-1;
		}
		double[] co = new double[coe.length];
		System.arraycopy(coe, 0, co, 0, coe.length-1);
		co[0] = coe[coe.length-1];
		return co;
	}

	/**
	 * Calculate score.
	 */
	public void CalculateScore() {
		Logistic Reg = new Logistic();
		try {
			Reg.buildClassifier(LS);
		} catch (Exception err) {
			err.printStackTrace();
		}
		// add by magic
		String s = Reg.toString();
		String[] ss = s.split("\n");
		
		String[] coeffi = new String[LS.numAttributes()+1-1];
		for (int i = 0; i < coeffi.length; i++) {
			coeffi[i] = ss[i + 3].trim();
			coeffi[i] = coeffi[i].split("\\s+")[1];
		}

		double[] co = stringTOdouble(coeffi);

		for( int i=0; i<co.length; i++)
		{
			co[i]*=-1;
//			The difinition of "1" and "0" in weka is different from that in R package;
			   //For example, for response(1,0,1,0,1), weka estimate of Intercept is -0.4055 (0.4055)
            //because it reversed 1 and 0 here.
            //In order to reverse it back, each co[] multiply by -1;
		}
		
/*
		co[0] = 1;
		co[1] = -5.29;
*/
		for( int j=0; j<Observe.getRowDimension(); j++)
		{
			double ex=0;
			for( int k=0; k<Observe.getColumnDimension(); k++)
			{
				ex+=Observe.get(j, k)*co[k];
			}
			ex+=co[co.length-1];
			double sc=Math.exp(ex)/(Math.exp(ex)+1);		
			Residual[j]=Response.get(j,0)-sc;
		}
	}

	/**
	 * Convert an array of strings to an array of double values. 
	 * @param s an array of strings, which are numerical as a matter of fact, respectively. 
	 * @return
	 */
	private static double[] stringTOdouble(String[] s) {
		if (s == null) {
			return null;
		}
//		for (int i=0; i < s.length; i++)
//		{
//			System.out.println(s[i]);
//		}
		double[] d = new double[s.length];
		for (int i = 0; i < d.length; i++) {
			d[i] = Double.parseDouble(s[i]);
		}
		return d;
	}

	public static void main(String[] args) {	
		double[][] vals = { { 0.8 }, { 0.3 }, { 0.2 }, { 0.6 }, { 0.7 } };
		Matrix A = new Matrix(vals);
		double[][] response = { { 1 }, { 0 }, { 1 }, { 0 }, { 1 } };
		Matrix Y = new Matrix(response);
		LogisticScore QS = new LogisticScore(vals, response);
		QS.CalculateScore();	
	}
}
