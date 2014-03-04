package gear.subcommands.metawatchdog.powercalculator;

import java.io.DataOutputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;

import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.FileUtil;
import gear.util.Logger;

public class DogPowerCommandImpl extends CommandImpl
{
	@Override
	public void execute(CommandArguments cmdArgs)
	{
		DogPowerCommandArguments dogpowerArgs = (DogPowerCommandArguments)cmdArgs;
		NormalDistribution nd = new NormalDistributionImpl();
		
		double alpha = dogpowerArgs.getAlpha();
		double beta = dogpowerArgs.getBeta();
		int tests = dogpowerArgs.getTests();
		
		double za = 0;
		double zb = 0;
		try
		{
			za = nd.inverseCumulativeProbability(1 - alpha / tests);
			zb = nd.inverseCumulativeProbability(1 - beta);
		}
		catch (MathException e)
		{
			Logger.handleException(e, "A math exception occurred when calculating inverse cumulative probability.");
		}
		
		double h2 = dogpowerArgs.getH2(); 

		double k = (zb + za/Math.sqrt(1-h2 * h2)) * (zb + za /Math.sqrt(1-h2 * h2)) * (1-h2 * h2) / (h2 * h2);
		Logger.printUserLog("alpha: " + alpha + " (z score=" + za + ")");
		Logger.printUserLog("beta: " + beta + " (z score=" + zb + ")");
		Logger.printUserLog("h2: " + h2);
		Logger.printUserLog("tests: " + tests);
		Logger.printUserLog("sample size to archive type I error rate at alphe = " + alpha + " and power = " + beta + " is " + k);

		PrintStream predictorFile = FileUtil.CreatePrintStream(dogpowerArgs.getOutRoot() + ".dogpower");
		predictorFile.println("alpha (type I error rate): " + alpha + "( z score =" + za + ")");
		predictorFile.println("beta (type II error rate): " + beta + "( z score =" + zb + ")");
		predictorFile.println("h2: " + h2);
		predictorFile.println("tests: " + tests);
		predictorFile.println("sample size to archive type I error rate at alphe = " + alpha + " and power = " + beta + " is " + k);

		DataOutputStream os = null;
		try
		{
			os = new DataOutputStream(new FileOutputStream(dogpowerArgs.getOutRoot() + ".encode"));
		}
		catch (FileNotFoundException e)
		{
			Logger.printUserError("Cannot create file '" + dogpowerArgs.getOutRoot()
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
			os.writeLong(dogpowerArgs.getSeed());
			os.writeDouble(k);
			os.close();
		}
		catch (IOException e)
		{
			Logger.handleException(e, "An I/O exception occurred when writing ecode.");
		}
	}
}
