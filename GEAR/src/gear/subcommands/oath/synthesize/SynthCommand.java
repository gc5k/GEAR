package gear.subcommands.oath.synthesize;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class SynthCommand extends Command
{

	@Override
	public String getName()
	{
		return "oath";
	}

	@Override
	public String getDescription()
	{
		return "Open GWAS algorithm (OATH).";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options)
	{
		options.addOption(OptionBuilder.withDescription(OPT_META_BATCH_DESC).withLongOpt(OPT_META_BATCH_LONG).hasArg().create(OPT_META_BATCH));
		options.addOption(OptionBuilder.withDescription(OPT_META_GZ_BATCH_DESC).withLongOpt(OPT_META_GZ_BATCH_LONG).hasArg().create(OPT_META_GZ_BATCH));
		options.addOption(OptionBuilder.withDescription(OPT_KEEP_ATGC_DESC).withLongOpt(OPT_KEEP_ATGC_LONG).create(OPT_KEEP_ATGC));
		options.addOption(OptionBuilder.withDescription(OPT_CM_DESC).hasArg().isRequired().create(OPT_CM));
		options.addOption(OptionBuilder.withDescription(OPT_CHR_DESC).hasArg().create(OPT_CHR));
		options.addOption(OptionBuilder.withDescription(OPT_N_DESC).hasArg().create(OPT_N));
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException
	{
		SynthCommandArguments synArgs = new SynthCommandArguments();

		if (cmdLine.hasOption(OPT_META_BATCH))
		{
			synArgs.setMetaBatch(cmdLine.getOptionValue(OPT_META_BATCH));
			synArgs.setGZ(false);
		}

		if (cmdLine.hasOption(OPT_META_GZ_BATCH))
		{
			synArgs.setMetaBatch(cmdLine.getOptionValue(OPT_META_GZ_BATCH));
			synArgs.setGZ(true);
		}

		synArgs.setCMFile(cmdLine.getOptionValue(OPT_CM));
		if (cmdLine.hasOption(OPT_CHR))
		{
			synArgs.setChr(cmdLine.getOptionValue(OPT_CHR));
		}
		
		if (cmdLine.hasOption(OPT_N))
		{
			synArgs.setN(cmdLine.getOptionValue(OPT_N));
		}

		return synArgs;
	}

	@Override
	protected CommandImpl createCommandImpl() 
	{
		SynthCommandImpl synImpl = new SynthCommandImpl();
		return synImpl;
	}

	private final static String OPT_META_BATCH = "mb";
	private final static String OPT_META_BATCH_LONG = "meta-batch";
	private final static String OPT_META_BATCH_DESC = "The summary statistic batch";

	private final static String OPT_META_GZ_BATCH = "mgb";
	private final static String OPT_META_GZ_BATCH_LONG = "meta-gz-batch";
	private final static String OPT_META_GZ_BATCH_DESC = "The summary statistic files in gz format";

	private final static String OPT_CHR = "chr";
	private final static String OPT_CHR_DESC = "Choose chromosome for analysis";

	private final static String OPT_KEEP_ATGC = "k";
	private final static String OPT_KEEP_ATGC_LONG = "keep-atgc";
	private final static String OPT_KEEP_ATGC_DESC = "Keep parlindormic loci.";

	private final static String OPT_CM = "cm";
	private final static String OPT_CM_DESC = "Correlation matrix";
	
	private final static String OPT_N = "n";
	private final static String OPT_N_DESC = "Sample size";

}
