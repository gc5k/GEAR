package gear.subcommands.he.assocpower;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;

import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.Logger;

public class HEAssocPowerCommandImpl extends CommandImpl
{
	@Override
	public void execute(CommandArguments cmdArgs)
	{
		try
		{
			HEAssocPowerCommandArguments heArgs = (HEAssocPowerCommandArguments)cmdArgs;
			NormalDistribution norm = new NormalDistributionImpl();
			double q1 = norm.inverseCumulativeProbability(1 - heArgs.getPower());
			double q2 = norm.inverseCumulativeProbability(1 - heArgs.getAlpha());
			double u = (q1 + q2) / heArgs.getHeritabilityVariance();
			long sampleSize = Math.round(2 * u * Math.sqrt(heArgs.getNumberOfMarkers()));
			Logger.printUserLog("Required sample size: " + sampleSize);
		}
		catch (MathException e)
		{
			Logger.handleException(e, "A math exception occured.");
		}
	}
}
