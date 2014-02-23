package gear.metawatchdog.powercalculator;

import java.io.DataOutputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;

import gear.CmdArgs;
import gear.util.FileUtil;
import gear.util.Logger;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.*;

public class MetaWatchdogPowerCalculator
{
	private int tests = 100;
	private double alpha = 0.05;
	private double beta = 0.5;
	private double h2 = 0.95;
	private long eseed;
	private NormalDistribution nd = new NormalDistributionImpl();
	
	public MetaWatchdogPowerCalculator()
	{
		alpha = CmdArgs.INSTANCE.dog_alpha;
		beta = CmdArgs.INSTANCE.dog_beta;
		tests = CmdArgs.INSTANCE.dog_tests;
		h2 = CmdArgs.INSTANCE.dog_h2;
		eseed = CmdArgs.INSTANCE.simuSeed;
		
		double za = 0;
		double zb = 0;
		try
		{
			za = nd.inverseCumulativeProbability(1-alpha/tests);
			zb = nd.inverseCumulativeProbability(1-beta);
		}
		catch (MathException e)
		{
			e.printStackTrace();
		}

		double k = (zb + za/Math.sqrt(1-h2 * h2)) * (zb + za /Math.sqrt(1-h2 * h2)) * (1-h2 * h2) / (h2 * h2);
		Logger.printUserLog("alpha: " + alpha + " (z score=" + za + ")");
		Logger.printUserLog("beta: " + beta + " (z score=" + zb + ")");
		Logger.printUserLog("h2: " + h2);
		Logger.printUserLog("tests: " + tests);
		Logger.printUserLog("sample size to archive type I error rate at alphe = " + alpha + " and power = " + beta + " is " + k);

		PrintStream predictorFile = FileUtil.CreatePrintStream(CmdArgs.INSTANCE.out + ".dogpower");
		predictorFile.println("alpha (type I error rate): " + alpha + "( z score =" + za + ")");
		predictorFile.println("beta (type II error rate): " + beta + "( z score =" + zb + ")");
		predictorFile.println("h2: " + h2);
		predictorFile.println("tests: " + tests);
		predictorFile.println("sample size to archive type I error rate at alphe = " + alpha + " and power = " + beta + " is " + k);

		DataOutputStream os = null;
		try
		{
			os = new DataOutputStream(new FileOutputStream(CmdArgs.INSTANCE.out + ".encode"));
		}
		catch (FileNotFoundException e)
		{
			Logger.printUserError("Cannot create file '" + CmdArgs.INSTANCE.out
					+ "'.");
			Logger.printUserError("Exception Message: " + e.getMessage());
			System.exit(1);
		}
		try
		{
			os.writeDouble(alpha);
			os.writeDouble(beta);
			os.writeInt(tests);
			os.writeDouble(h2);
			os.writeLong(eseed);
			os.writeDouble(k);
			os.close();
		}
		catch (IOException e)
		{
			Logger.handleException(e, "An I/O exception occurred when writing ecode.");
		}
	}
}
