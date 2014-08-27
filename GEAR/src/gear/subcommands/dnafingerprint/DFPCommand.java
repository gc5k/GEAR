package gear.subcommands.dnafingerprint;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class DFPCommand extends Command
{
	@Override
	public String getName()
	{
		return "dfp";
	}

	@Override
	public String getDescription()
	{
		return "Genotype similarity analysis";
	}

	@Override
	public void prepareOptions(Options options)
	{
		options.addOption(OptionBuilder.withDescription(OPT_BFILE2_DESC).withLongOpt(OPT_BFILE2_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_NUM_MARKER_DESC).withLongOpt(OPT_NUM_MARKER_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_LOW_CUTOFF_DESC).withLongOpt(OPT_LOW_CUTOFF_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_HIGH_CUTOFF_DESC).withLongOpt(OPT_HIGH_CUTOFF_LONG).hasArg().create());
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine)
			throws CommandArgumentException
	{
		DFPCommandArguments cmdArgs = new DFPCommandArguments();
		if (cmdLine.hasOption(OPT_BFILE2_LONG))
		{
			cmdArgs.setBFile2(cmdLine.getOptionValue(OPT_BFILE2_LONG));
		}
		
		cmdArgs.setLowCutoff(parseDoubleOptionValueInRange(cmdLine, OPT_LOW_CUTOFF_LONG, OPT_LOW_CUTOFF_DEFAULT, -5, 1));
		cmdArgs.setHighCutoff(parseDoubleOptionValueInRange(cmdLine, OPT_HIGH_CUTOFF_LONG, OPT_HIGH_CUTOFF_DEFAULT, 0, 1));
		cmdArgs.setNumMarker(parseLongOptionValueInRange(cmdLine, OPT_NUM_MARKER_LONG, OPT_NUM_MARKER_DEFAULT, 0, Long.MAX_VALUE));
		return null;
	}

	@Override
	protected CommandImpl createCommandImpl()
	{
		
		return new DFPCommandImpl();
	}

	private static final String OPT_BFILE2_LONG = "bfile2";
	private static final String OPT_BFILE2_DESC = "The second bfile.";
	
	private static final String OPT_NUM_MARKER_LONG = "num-marker";
	private static final String OPT_NUM_MARKER_DEFAULT = "-1";
	private static final String OPT_NUM_MARKER_DESC = "The number of markers for comparison.";
	
	private static final String OPT_LOW_CUTOFF_LONG = "low-cutoff";
	private static final String OPT_LOW_CUTOFF_DEFAULT = "-1";
	private static final String OPT_LOW_CUTOFF_DESC = "The lower bound cutoff for similarity. Default is -1.0";

	private static final String OPT_HIGH_CUTOFF_LONG = "high-cutoff";
	private static final String OPT_HIGH_CUTOFF_DEFAULT = "1";
	private static final String OPT_HIGH_CUTOFF_DESC = "The higher bound cutoff for similarity. Default is 1.0";
}
