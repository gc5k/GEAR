package gear.subcommands.metawatchdog.powercalculator;

import java.io.DataOutputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ChiSquaredDistributionImpl;
import org.apache.commons.math.distribution.NormalDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;

import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.Logger;

public class DogPowerCommandImpl extends CommandImpl
{
	private DogPowerCommandArguments dogpowerArgs;
	private int method = 0; //0 for chisq; 1 for regression;
	@Override
	public void execute(CommandArguments cmdArgs)
	{
		dogpowerArgs = (DogPowerCommandArguments) cmdArgs;

		if(dogpowerArgs.getChisqFlag())
		{
			method = 0;
			chisq();
		}
		else if (dogpowerArgs.getRegressionFlag())
		{
			method = 1;
			regression();
		}
		else
		{
			Logger.printUserLog("No method has been specified. GEAR quitted.");
			System.exit(0);
		}
		
		Logger.printUserLog("Parameters have been saved into '" + dogpowerArgs.getOutRoot() + ".encode'.");
	}

	private void chisq()
	{
		ChiSquaredDistributionImpl chiDis = new ChiSquaredDistributionImpl(1);

		double alpha = dogpowerArgs.getAlpha();
		int tests = dogpowerArgs.getTests();
		double q = dogpowerArgs.getQValue();
		double logP = -1 * Math.log10(alpha/tests);

		double logChisqP = 0;
		try
		{
			logChisqP = -1 * Math.log10(chiDis.cumulativeProbability(q));
		}
		catch (MathException e)
		{
			Logger.handleException(e, "error in getting p value for chisq distribution.");
		}
		double dif = logChisqP - logP;
		double k = 1;

		while (dif < 0)
		{
			k = k + 1;
			chiDis = new ChiSquaredDistributionImpl(k);
			try
			{
				logChisqP = -1 * Math.log10(chiDis.cumulativeProbability(q));
			}
			catch (MathException e)
			{
				Logger.handleException(e, "error in getting p value for chisq distribution.");
			}
			dif = logChisqP - logP;
		}
		
		Logger.printUserLog("Alpha: " + alpha);
		Logger.printUserLog("Q value: " + q);
		Logger.printUserLog("Tests: " + tests);
		Logger.printUserLog("Method: chisq");
		Logger.printUserLog("Sample size to archive type I error rate at alphe = " + alpha + " is " + k);

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
			os.writeDouble(k);
			os.writeLong(dogpowerArgs.getSeed());
			os.writeDouble(dogpowerArgs.getAlpha());
			os.writeInt(dogpowerArgs.getTests());

			os.writeDouble(dogpowerArgs.getBeta());
			os.writeDouble(dogpowerArgs.getB());
			os.writeDouble(dogpowerArgs.getQValue());
			os.writeInt(method);

			os.close();
		}
		catch (IOException e)
		{
			Logger.handleException(e, "An I/O exception occurred when writing ecode.");
		}
	}

	private void regression()
	{
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

		double b = dogpowerArgs.getB();

		double k = (zb + za/Math.sqrt(1-b * b)) * (zb + za /Math.sqrt(1-b * b)) * (1-b * b) / (b * b);
		Logger.printUserLog("Alpha: " + alpha);
		Logger.printUserLog("Beta: " + beta);
		Logger.printUserLog("Regression coefficient: " + b);
		Logger.printUserLog("Tests: " + tests);
		Logger.printUserLog("Method: regression");
		Logger.printUserLog("Sample size to archive type I error rate at alphe = " + alpha + " and power = " + beta + " is " + k);

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
			os.writeDouble(k);
			os.writeLong(dogpowerArgs.getSeed());
			os.writeDouble(dogpowerArgs.getAlpha());
			os.writeInt(dogpowerArgs.getTests());

			os.writeDouble(dogpowerArgs.getBeta());
			os.writeDouble(dogpowerArgs.getB());
			os.writeDouble(dogpowerArgs.getQValue());
			os.writeInt(method);

			os.close();
		}
		catch (IOException e)
		{
			Logger.handleException(e, "An I/O exception occurred when writing ecode.");
		}
	}
}
