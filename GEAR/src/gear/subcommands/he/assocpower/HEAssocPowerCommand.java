package gear.subcommands.he.assocpower;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class HEAssocPowerCommand extends Command
{
	@Override
	public String getName()
	{
		return "he-assoc-power";
	}

	@Override
	public String getDescription()
	{
		return "Haseman-Elston Association Power Calculator";
	}
	
	@Override
	public String getLongDescription()
	{
		return "Calculate how many samples are needed to reach a given power of Haseman-Elston association.";
	}

	@SuppressWarnings("static-access")
	@Override
	protected void prepareOptions(Options options)
	{
		options.addOption(OptionBuilder.withDescription(OPT_NUMBER_OF_MARKERS_DESC).withLongOpt(OPT_NUMBER_OF_MARKERS_LONG).hasArg().isRequired().create(OPT_NUMBER_OF_MARKERS));
		options.addOption(OptionBuilder.withDescription(OPT_HERITABILITY_VARIANCE_DESC).withLongOpt(OPT_HERITABILITY_VARIANCE_LONG).hasArg().isRequired().create(OPT_HERITABILITY_VARIANCE));
		options.addOption(OptionBuilder.withDescription(OPT_ALPHA_DESC).withLongOpt(OPT_ALPHA_LONG).hasArg().create(OPT_ALPHA));
		options.addOption(OptionBuilder.withDescription(OPT_POWER_DESC).withLongOpt(OPT_POWER_LONG).hasArg().isRequired().create(OPT_POWER));
	}

	@Override
	protected CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException
	{
		HEAssocPowerCommandArguments cmdArgs = new HEAssocPowerCommandArguments();
		
		try
		{
			cmdArgs.setNumberOfMarkers(Integer.valueOf(cmdLine.getOptionValue(OPT_NUMBER_OF_MARKERS)));
		}
		catch (NumberFormatException e)
		{
			throw new CommandArgumentException("'" + cmdLine.getOptionValue(OPT_NUMBER_OF_MARKERS) + "' is not a valid integer.");
		}
		
		String optValue = "";
		try
		{
			optValue = cmdLine.getOptionValue(OPT_HERITABILITY_VARIANCE);
			cmdArgs.setHeritabilityVariance(Float.valueOf(optValue));
			
			optValue = cmdLine.getOptionValue(OPT_ALPHA, OPT_ALPHA_DEFAULT);
			float alpha = Float.valueOf(optValue);
			if (alpha < 0 || alpha > 1.0)
			{
				throw new CommandArgumentException("Alpha must lie between 0 and 1");
			}
			cmdArgs.setAlpha(alpha);
			
			optValue = cmdLine.getOptionValue(OPT_POWER);
			float power = Float.valueOf(optValue);
			if (power < 0 || power > 1.0)
			{
				throw new CommandArgumentException("Power must lie between 0 and 1");
			}
			cmdArgs.setPower(power);
		}
		catch (NumberFormatException e)
		{
			throw new CommandArgumentException("'" + optValue + "' is not a valid floating point number.");
		}
		
		return cmdArgs;
	}

	@Override
	protected CommandImpl createCommandImpl()
	{
		return new HEAssocPowerCommandImpl();
	}
	
	private static final char OPT_NUMBER_OF_MARKERS = 'm';
	private static final String OPT_NUMBER_OF_MARKERS_LONG = "markers";
	private static final String OPT_NUMBER_OF_MARKERS_DESC = "effective number of markers";
	
	private static final char OPT_HERITABILITY_VARIANCE = 'h';
	private static final String OPT_HERITABILITY_VARIANCE_LONG = "heritability2";
	private static final String OPT_HERITABILITY_VARIANCE_DESC = "sampling variance of heritability";
	
	private static final String OPT_ALPHA_DEFAULT = "0.05";
	private static final char OPT_ALPHA = 'a';
	private static final String OPT_ALPHA_LONG = "alpha";
	private static final String OPT_ALPHA_DESC = "type I error rate, default to " + OPT_ALPHA_DEFAULT;
	
	private static final char OPT_POWER = 'p';
	private static final String OPT_POWER_LONG = "power";
	private static final String OPT_POWER_DESC = "statistical power";
}
