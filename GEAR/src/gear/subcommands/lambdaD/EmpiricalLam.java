package gear.subcommands.lambdaD;

import gear.util.Logger;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ChiSquaredDistributionImpl;
import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;

public class EmpiricalLam
{
	public EmpiricalLam(double[] DesStat, double size1, double size2)
	{
		this.DesStat = DesStat;
		this.size1 = size1;
		this.size2 = size2;
		CalculateLamD();
		QT();
	}

	public EmpiricalLam(double[] DesStat, double cs1, double ctrl1, double cs2, double ctrl2)
	{
		this.DesStat = DesStat;
		this.cs1 = cs1;
		this.ctrl1 = ctrl1;
		this.cs2 = cs2;
		this.ctrl2 = ctrl2;
		CalculateLamD();
		CC();
	}

	private void CalculateLamD()
	{
		ChiSquaredDistributionImpl chiDis = new ChiSquaredDistributionImpl(1);
		for (int i = (int) (DesStat.length * 0.05); i < (int) (DesStat.length * 0.95); i++)
		{ //only use the points in 5%~95% interval
			try
			{
				double chi0 = chiDis.inverseCumulativeProbability((i+1)/(DesStat.length+0.05));
				Lam.addValue(DesStat[i]/chi0);
			}
			catch (MathException e)
			{
				e.printStackTrace();
			}
		}
	}

	private void QT()
	{
		double k = 2/(Math.sqrt(size1 / size2) + Math.sqrt(size1 / size2));		
		EmpLamMean = Lam.getMean();
		EmpLamSD = Lam.getStandardDeviation();

		EmpRhoMean = (1-EmpLamMean)/k;
		EmpRhoSD = EmpLamSD/k;

		EmpOSMean = EmpRhoMean * Math.sqrt(size1 * size2);
		EmpOSSD = EmpLamSD/k * Math.sqrt(size1 * size2);
	}

	private void CC()
	{
		double k = 2/(Math.sqrt((cs1 + ctrl1) / (cs2 + ctrl2)) + Math.sqrt( (cs2 + ctrl2) / (cs1 + ctrl1)) );
		double R1 = cs1/ctrl1;
		double R2 = cs2/ctrl2;

		EmpLamMean = Lam.getMean();
		EmpLamSD = Lam.getStandardDeviation();
		
		EmpRhoMean = (1-EmpLamMean)/k;
		EmpRhoSD = EmpLamSD/k;

		EmpOSMean = (1-EmpLamMean)/k * Math.sqrt( (cs1 + ctrl1) * (cs2 + ctrl2) );
		EmpOSSD = EmpLamSD/k * Math.sqrt( (cs1 + ctrl1) * (cs2 + ctrl2) );
		
		EmpOSCtrlMean = (1-EmpLamMean) /k * Math.sqrt( (cs1 + ctrl1) * (cs2 + ctrl2) ) / Math.sqrt(R1 * R2);
		EmpOSCtrlSD = EmpLamSD/k * Math.sqrt( (cs1 + ctrl1) * (cs2 + ctrl2) ) / Math.sqrt(R1 * R2);
		
		EmpOSCsMean = (1-EmpLamMean) /k * Math.sqrt( (cs1 + ctrl1) * (cs2 + ctrl2) ) * Math.sqrt(R1 * R2);
		EmpOSCsSD = EmpLamSD/k * Math.sqrt( (cs1 + ctrl1) * (cs2 + ctrl2) ) * Math.sqrt(R1 * R2);	
	}

	protected void PrintCC()
	{
		Logger.printUserLog("Empirical Distribution for lambdaD: ");
		Logger.printUserLog("Mean: " + EmpLamMean);
		Logger.printUserLog("Standard deviation: " + EmpLamSD + "\n");

		Logger.printUserLog("Empirical Distribution for correlation:");
		Logger.printUserLog("Mean: " + EmpRhoMean);
		Logger.printUserLog("Standard deviation: " + EmpRhoSD + "\n");

		Logger.printUserLog("Empirical Overlapping Samples: ");
		Logger.printUserLog("Overlapping samples: " + EmpOSMean + " (" + EmpOSSD + ")");
		Logger.printUserLog("Overlapping control samples: " + EmpOSCtrlMean + " (" + EmpOSCtrlSD + ")");
		Logger.printUserLog("Overlapping case samples: " + EmpOSCsMean + " (" + EmpOSCsSD + ")\n");
	}

	protected void PrintQT()
	{
		Logger.printUserLog("Empirical Distribution for lambdaD:");
		Logger.printUserLog("Mean: " + EmpLamMean);
		Logger.printUserLog("Standard deviation: " + EmpLamSD);

		Logger.printUserLog("Empirical Overlapping Samples: ");
		Logger.printUserLog("Overlapping samples: " + EmpOSMean + " (" + EmpOSSD + ")");
	}
	
	private double[] DesStat;
	private DescriptiveStatistics Lam = new DescriptiveStatistics();

	private double EmpLamMean = 0;
	private double EmpLamSD = 0;

	private double EmpRhoMean = 0;
	private double EmpRhoSD = 0;
	
	private double EmpOSMean = 0;
	private double EmpOSSD = 0;

	private double EmpOSCtrlMean = 0;
	private double EmpOSCtrlSD = 0;

	private double EmpOSCsMean = 0;
	private double EmpOSCsSD = 0;

	private double size1 = 0;
	private double size2 = 0;

	private double cs1 = 0;
	private double ctrl1 = 0;
	private double cs2 = 0;
	private double ctrl2 = 0;
}
