package gear.subcommands.oath.oathx;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class OATHXCommand extends Command 
{

	@Override
	public String getName() 
	{
		return "oathX";
	}

	@Override
	public String getDescription() 
	{
		return "Exhaustive OATH";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options)
	{
		options.addOption(OptionBuilder.withDescription(OPT_NSS_BATCH_DESC).withLongOpt(OPT_NSS_BATCH_LONG).hasArg().isRequired().create(OPT_NSS_BATCH));
		options.addOption(OptionBuilder.withDescription(OPT_CM_DESC).hasArg().isRequired().create(OPT_CM));
		options.addOption(OptionBuilder.withDescription(OPT_N_DESC).hasArg().isRequired().create(OPT_N));
		options.addOption(OptionBuilder.withDescription(OPT_VERBOSE_LONG_DESC).withLongOpt(OPT_VERBOSE_LONG).create(OPT_VERBOSE));
		
		options.addOption(OptionBuilder.withDescription(OPT_CHR_DESC).withLongOpt(OPT_CHR_LONG).hasArgs().create());
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException 
	{
		OATHXCommandArguments obArgs = new OATHXCommandArguments();

		obArgs.setNSSBatch(cmdLine.getOptionValue(OPT_NSS_BATCH));
		obArgs.setCMFile(cmdLine.getOptionValue(OPT_CM));
		obArgs.setN(cmdLine.getOptionValue(OPT_N));

		if (cmdLine.hasOption(OPT_CHR_LONG))
		{
			obArgs.setChr(cmdLine.getOptionValues(OPT_CHR_LONG));
		}
		if (cmdLine.hasOption(OPT_VERBOSE_LONG))
		{
			obArgs.setVerbose();
		}

		return obArgs;
	}

	@Override
	protected CommandImpl createCommandImpl() 
	{
		return new OATHXCommandImpl();
	}

	private final static String OPT_NSS_BATCH = "nb";
	private final static String OPT_NSS_BATCH_LONG = "nss-batch";
	private final static String OPT_NSS_BATCH_DESC = "The summary statistic batch";

	private final static String OPT_CM = "cm";
	private final static String OPT_CM_DESC = "Correlation matrix";

	private static final String OPT_N = "n";
	private static final String OPT_N_DESC = "Specify the sample size";

	private final static String OPT_VERBOSE = "v";
	private final static String OPT_VERBOSE_LONG  = "verbose";
	private final static String OPT_VERBOSE_LONG_DESC  = "verbose";

}
