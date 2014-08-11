package gear.subcommands.lambdaD;

import gear.util.Logger;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ChiSquaredDistributionImpl;
import org.apache.commons.math.distribution.NormalDistributionImpl;
import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;

public class XTest
{

	public XTest(double[] DesStat, double n1, double n2)
	{
		this.DesStat = DesStat;
		this.n1 = n1;
		this.n2 = n2;
		Me = DesStat.length;
		CalX();
		QT();
	}

	public XTest(double[] DesStat, double cs1, double ctrl1, double cs2, double ctrl2)
	{
		this.DesStat = DesStat;
		this.cs1 = cs1;
		this.ctrl1 = ctrl1;
		this.cs2 = cs2;
		this.ctrl2 = ctrl2;
		Me = DesStat.length;
		CalX();
//		CC();
	}

	private void CalX()
	{
		for(int i = 0; i < DesStat.length; i++)
		{
			XVec.addValue(DesStat[i]);
		}
		ChiSquaredDistributionImpl chiDis = new ChiSquaredDistributionImpl(Me);

		try
		{
			Chi0025 = chiDis.inverseCumulativeProbability(0.025);
			Chi005 = chiDis.inverseCumulativeProbability(0.05);
			Chi095 = chiDis.inverseCumulativeProbability(0.95);
			Chi0975 = chiDis.inverseCumulativeProbability(0.975);
		}
		catch (MathException e)
		{
			Logger.handleException(e, "error in getting pvalue for Chi-sq in lambda.");
		}

		try
		{
			double z= (XVec.getSum() - Me)/Math.sqrt(2*Me);
			pX = 2* (1-nDis.cumulativeProbability(Math.abs(z)));
		}
		catch (MathException e)
		{
			Logger.handleException(e, "error in getting pvalue.");
		}
		
	}

	private void QT()
	{
		rho = (1-XVec.getSum()/Me) * (Math.sqrt(n1/n2) + Math.sqrt(n2/n1))/2;
		n12 = (1-XVec.getSum()/Me) * (n1+n2)/2;
		sigma_n12 = (n1+n2)/Math.sqrt(2*Me);

		try
		{
			pN12 = 2 * (1-nDis.cumulativeProbability(Math.abs(n12/sigma_n12)));
		}
		catch (MathException e)
		{
			Logger.handleException(e, "error in getting pvalue.");
		}
	}

	protected void PrintQT()
	{
		Logger.printUserLog("X statistic: " + XVec.getSum()/Me);
		Logger.printUserLog("p-value for X statistic: " + pX);
		Logger.printUserLog("95% CI for overlapping samples: "  + "(" + (-1.96*Math.sqrt(2/Me)) + "," + 1.96*Math.sqrt(2/Me) + ")");		
		Logger.printUserLog("Overlapping Samples: " + n12 + "\n");
		Logger.printUserLog("p-value for overlapping sample: " + pN12 + "\n");
		Logger.printUserLog("95% CI for overlapping samples: "  + "(" + -1.96*sigma_n12 + "," + 1.96*sigma_n12 + ")");
	}

	protected double getX()
	{
		return XVec.getSum()/Me;
	}

	protected double getpX()
	{
		return pX;
	}

	protected double getN12()
	{
		return n12;
	}

	protected double getrho()
	{
		return rho;
	}
	
	private ChiSquaredDistributionImpl chiDis = null;
	private NormalDistributionImpl nDis = new NormalDistributionImpl();

	private double[] DesStat;
	private DescriptiveStatistics XVec = new DescriptiveStatistics();

	private double Me = 0;
	private double pX = 0;
	private double rho = 0;
	
	private DescriptiveStatistics Lam = new DescriptiveStatistics();

	private double n1 = 0;
	private double n2 = 0;

	private double lambda = 0;

	private double cs1 = 0;
	private double ctrl1 = 0;
	private double cs2 = 0;
	private double ctrl2 = 0;

	private double Chi0025 = 0;
	private double Chi005 = 0;
	private double Chi0975 = 0;
	private double Chi095 = 0;

	private double n12 = 0;
	private double sigma_n12 = 0;
	private double pN12 = 0;

	private double n12L = 0;
	private double n12H = 0;
}
