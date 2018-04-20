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
	
	public DFPCommand()
	{
		addAlias("dfp");
	}

	@Override
	public String getName()
	{
		return "dnafinger";
	}

	@Override
	public String getDescription()
	{
		return "Genotype similarity analysis";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options)
	{
		options.addOption(OptionBuilder.withDescription(OPT_BFILE_DESC).withLongOpt(OPT_BFILE_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_BFILE2_DESC).withLongOpt(OPT_BFILE2_LONG).hasArg().create());

		options.addOption(OptionBuilder.withDescription(OPT_NUM_MARKER_DESC).withLongOpt(OPT_NUM_MARKER_LONG).hasArg().create(OPT_NUM_MARKER));
		options.addOption(OptionBuilder.withDescription(OPT_LOW_CUTOFF_DESC).withLongOpt(OPT_LOW_CUTOFF_LONG).hasArg().create(OPT_LOW_CUTOFF));
		options.addOption(OptionBuilder.withDescription(OPT_HIGH_CUTOFF_DESC).withLongOpt(OPT_HIGH_CUTOFF_LONG).hasArg().create(OPT_HIGH_CUTOFF));
		options.addOption(OptionBuilder.withDescription(OPT_SEED_DESC).withLongOpt(OPT_SEED_LONG).hasArg().create());
		
	    options.addOption(OptionBuilder.withDescription(OPT_KEEP_DESC).withLongOpt(OPT_KEEP_LONG).hasArg().create());
	    options.addOption(OptionBuilder.withDescription(OPT_REMOVE_DESC).withLongOpt(OPT_REMOVE_LONG).hasArg().create());

	    options.addOption(OptionBuilder.withDescription(OPT_EXTRACT_DESC).withLongOpt(OPT_EXTRACT_LONG).hasArg().create());
	    options.addOption(OptionBuilder.withDescription(OPT_EXCLUDE_DESC).withLongOpt(OPT_EXCLUDE_LONG).hasArg().create());

	    options.addOption(OptionBuilder.withDescription(OPT_CHR_DESC).withLongOpt(OPT_CHR_LONG).hasArgs().create());
	    options.addOption(OptionBuilder.withDescription(OPT_NOT_CHR_DESC).withLongOpt(OPT_NOT_CHR_LONG).hasArgs().create());

	}

	@Override
	public CommandArguments parse(CommandLine cmdLine)
			throws CommandArgumentException
	{
		DFPCommandArguments dfpArgs = new DFPCommandArguments();

	    parseFileArguments((CommandArguments) dfpArgs, cmdLine);
	    parseSampleFilterArguments((CommandArguments) dfpArgs, cmdLine);
	    parseSNPFilterChromosomeArguments((CommandArguments) dfpArgs, cmdLine);

		if (cmdLine.hasOption(OPT_BFILE2_LONG))
		{
			dfpArgs.setBFile2(cmdLine.getOptionValue(OPT_BFILE2_LONG));
		}

		dfpArgs.setLowCutoff(parseDoubleOptionValueInRange(cmdLine, OPT_LOW_CUTOFF_LONG, OPT_LOW_CUTOFF_DEFAULT, -5, 1));
		dfpArgs.setHighCutoff(parseDoubleOptionValueInRange(cmdLine, OPT_HIGH_CUTOFF_LONG, OPT_HIGH_CUTOFF_DEFAULT, 0, 1));
		dfpArgs.setNumMarker(parseLongOptionValueInRange(cmdLine, OPT_NUM_MARKER_LONG, OPT_NUM_MARKER_DEFAULT, -1, Long.MAX_VALUE));
		dfpArgs.setSeed(parseLongOptionValueInRange(cmdLine, OPT_SEED_LONG, OPT_SEED_DEFAULT, 0, Long.MAX_VALUE));
		return dfpArgs;
	}

	@Override
	protected CommandImpl createCommandImpl()
	{
		
		return new DFPCommandImpl();
	}

	private static final String OPT_BFILE2_LONG = "bfile2";
	private static final String OPT_BFILE2_DESC = "The second bfile.";

	private static final String OPT_NUM_MARKER = "m";
	private static final String OPT_NUM_MARKER_LONG = "num-marker";
	private static final String OPT_NUM_MARKER_DEFAULT = "100";
	private static final String OPT_NUM_MARKER_DESC = "The number of markers for comparison.";
	
	private static final String OPT_LOW_CUTOFF = "l";
	private static final String OPT_LOW_CUTOFF_LONG = "low-cutoff";
	private static final String OPT_LOW_CUTOFF_DEFAULT = "-1";
	private static final String OPT_LOW_CUTOFF_DESC = "The lower bound cutoff for similarity. Default is -1.0";

	private static final String OPT_HIGH_CUTOFF = "h";
	private static final String OPT_HIGH_CUTOFF_LONG = "high-cutoff";
	private static final String OPT_HIGH_CUTOFF_DEFAULT = "1";
	private static final String OPT_HIGH_CUTOFF_DESC = "The higher bound cutoff for similarity. Default is 1.0";
}
