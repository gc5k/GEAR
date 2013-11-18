package gear.metawatchdog.powercalculator;

import gear.CmdArgs;
import gear.util.Logger;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.*;

public class MetaWatchdogPowerCalculator
{
	private int tests = 100;
	private double alpha = 0.05;
	private double beta = 0.9;
	private double h2 = 0.95;
	private NormalDistribution nd = new NormalDistributionImpl();
	
	public MetaWatchdogPowerCalculator()
	{
		alpha = CmdArgs.INSTANCE.dog_alpha;
		beta = CmdArgs.INSTANCE.dog_beta;
		tests = CmdArgs.INSTANCE.dog_tests;
		h2 = CmdArgs.INSTANCE.dog_h2;
		
		double za = 0;
		double zb = 0;
		try
		{
			za = nd.inverseCumulativeProbability(alpha/tests);
			zb = nd.inverseCumulativeProbability(beta);
		}
		catch (MathException e)
		{
			e.printStackTrace();
		}

		double k = (zb - za/Math.sqrt(1-h2 * h2)) * (zb - za /Math.sqrt(1-h2 * h2)) * (1-h2 * h2) / (h2 * h2);
		Logger.printUserLog("alpha: " + alpha + "(" + za + ")");
		Logger.printUserLog("bete: " + beta + "(" + zb + ")");
		Logger.printUserLog("h2: " + h2);
		Logger.printUserLog("tests: " + tests);
		Logger.printUserLog("sample size to archive type I error rate at alphe = " + alpha + " and power = " + beta + " is " + k);
	}
}
