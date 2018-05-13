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
		options.addOption(OptionBuilder.withDescription(OPT_NSS_BATCH_DESC).withLongOpt(OPT_NSS_BATCH_LONG).hasArg().create(OPT_NSS_BATCH));
		options.addOption(OptionBuilder.withDescription(OPT_NSS_GZ_BATCH_DESC).withLongOpt(OPT_NSS_GZ_BATCH_LONG).hasArg().create(OPT_NSS_GZ_BATCH));
//		options.addOption(OptionBuilder.withDescription(OPT_KEEP_ATGC_DESC).withLongOpt(OPT_KEEP_ATGC_LONG).create(OPT_KEEP_ATGC));
		options.addOption(OptionBuilder.withDescription(OPT_CM_DESC).hasArg().isRequired().create(OPT_CM));
		options.addOption(OptionBuilder.withDescription(OPT_N_DESC).hasArg().create(OPT_N));
		options.addOption(OptionBuilder.withDescription(OPT_KEEP_BATCH_DESC).withLongOpt(OPT_KEEP_BATCH_LONG).hasArgs().create());
		options.addOption(OptionBuilder.withDescription(OPT_VERBOSE_LONG_DESC).withLongOpt(OPT_VERBOSE_LONG).create(OPT_VERBOSE));
		
		options.addOption(OptionBuilder.withDescription(OPT_CHR_DESC).withLongOpt(OPT_CHR_LONG).hasArgs().create());
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException
	{
		SynthCommandArguments synArgs = new SynthCommandArguments();

		if (cmdLine.hasOption(OPT_NSS_BATCH))
		{
			synArgs.setNSSBatch(cmdLine.getOptionValue(OPT_NSS_BATCH));
			synArgs.setGZ(false);
		}

		if (cmdLine.hasOption(OPT_NSS_GZ_BATCH))
		{
			synArgs.setNSSBatch(cmdLine.getOptionValue(OPT_NSS_GZ_BATCH));
			synArgs.setGZ(true);
		}

		synArgs.setCMFile(cmdLine.getOptionValue(OPT_CM));
		if (cmdLine.hasOption(OPT_CHR_LONG))
		{
			synArgs.setChr(cmdLine.getOptionValues(OPT_CHR_LONG));
		}
		
		if (cmdLine.hasOption(OPT_N))
		{
			synArgs.setN(cmdLine.getOptionValue(OPT_N));
		}

		if (cmdLine.hasOption(OPT_KEEP_BATCH_LONG))
		{
			synArgs.setKeepBatch(cmdLine.getOptionValues(OPT_KEEP_BATCH_LONG));
		}
		
		if (cmdLine.hasOption(OPT_VERBOSE_LONG))
		{
			synArgs.setVerbose();
		}
		return synArgs;
	}

	@Override
	protected CommandImpl createCommandImpl() 
	{
		SynthCommandImpl synImpl = new SynthCommandImpl();
		return synImpl;
	}

	private final static String OPT_NSS_BATCH = "nb";
	private final static String OPT_NSS_BATCH_LONG = "nss-batch";
	private final static String OPT_NSS_BATCH_DESC = "The summary statistic batch";

	private final static String OPT_NSS_GZ_BATCH = "mgb";
	private final static String OPT_NSS_GZ_BATCH_LONG = "nss-gz-batch";
	private final static String OPT_NSS_GZ_BATCH_DESC = "The summary statistic files in gz format";

//	private final static String OPT_KEEP_ATGC = "k";
//	private final static String OPT_KEEP_ATGC_LONG = "keep-atgc";
//	private final static String OPT_KEEP_ATGC_DESC = "Keep parlindormic loci.";

	private final static String OPT_CM = "cm";
	private final static String OPT_CM_DESC = "Correlation matrix";
	
	private final static String OPT_N = "n";
	private final static String OPT_N_DESC = "Sample size";

	private final static String OPT_KEEP_BATCH_LONG = "keep-nss";
	private final static String OPT_KEEP_BATCH_DESC = "Specify the indexes of nss files.";

	private final static String OPT_VERBOSE = "v";
	private final static String OPT_VERBOSE_LONG  = "verbose";
	private final static String OPT_VERBOSE_LONG_DESC  = "verbose";
}
