package gear.subcommands.arch;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class ArchCommand extends Command
{

	@Override
	public String getName()
	{
		return "arch";
	}

	@Override
	public String getDescription()
	{
		return "Inference of genetic architecture";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options)
	{
//		options.addOption(OptionBuilder.withDescription(OPT_LD_R_LONG_DESC).withLongOpt(OPT_LD_R_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_BFILE_DESC).withLongOpt(OPT_BFILE_LONG).hasArg().isRequired().create());
		options.addOption(OptionBuilder.withDescription(OPT_QT_META_BATCH_LONG_DESC).withLongOpt(OPT_QT_META_BATCH_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_CC_META_BATCH_LONG_DESC).withLongOpt(OPT_CC_META_BATCH_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_KEY_DESC).hasArgs().create(OPT_KEY));
		options.addOption(OptionBuilder.withDescription(OPT_WINDOW_DESC).hasArgs().create(OPT_WINDOW));
		options.addOption(OptionBuilder.withDescription(OPT_EXTRACT_DESC).hasArgs().create(OPT_EXTRACT));
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException
	{
		ArchCommandArguments archArgs = new ArchCommandArguments();
		parseFileArguments(archArgs, cmdLine);

//		if (cmdLine.hasOption(OPT_LD_R_LONG))
//		{
//			archArgs.setLDFile(cmdLine.getOptionValue(OPT_LD_R_LONG));
//		}

		if (cmdLine.hasOption(OPT_WINDOW))
		{
			archArgs.setWindow(cmdLine.getOptionValue(OPT_WINDOW));
		}

		if (cmdLine.hasOption(OPT_QT_META_BATCH_LONG))
		{
			archArgs.setQtMetaBatch(cmdLine.getOptionValue(OPT_QT_META_BATCH_LONG));
		}

		if (cmdLine.hasOption(OPT_CC_META_BATCH_LONG))
		{
			archArgs.setCCMetaBatch(cmdLine.getOptionValue(OPT_CC_META_BATCH_LONG));
		}

		if (cmdLine.hasOption(OPT_EXTRACT))
		{
			archArgs.setExtract(cmdLine.getOptionValue(OPT_EXTRACT));
		}

		if (cmdLine.hasOption(OPT_KEY))
		{
			archArgs.setKey(cmdLine.getOptionValues(OPT_KEY));
		}

		return archArgs;
	}

	private void parseFileArguments(ArchCommandArguments archArgs, CommandLine cmdLine) throws CommandArgumentException
	{
		String bfile = cmdLine.getOptionValue("bfile");

		if (bfile == null)
		{
			throw new CommandArgumentException("No genotypes are provided. Either --bfile or --file must be set.");
		}

		archArgs.setBFile(bfile);
	}

	@Override
	protected CommandImpl createCommandImpl()
	{
		return new ArchCommandImpl();
	}
	
	private final static String OPT_WINDOW = "window";
	private final static String OPT_WINDOW_DESC = "SNP pairs having distance smaller than the values (mb) specified will have ld calculation";

//	private final static String OPT_LD_R_LONG = "ld-score";
//	private final static String OPT_LD_R_LONG_DESC = "Read ld from the file specified";

	private final static String OPT_QT_META_BATCH_LONG = "qt-meta-batch";
	private final static String OPT_QT_META_BATCH_LONG_DESC = "Read summary statistics for quantitative traits from the file specified";

	private final static String OPT_CC_META_BATCH_LONG = "cc-meta-batch";
	private final static String OPT_CC_META_BATCH_LONG_DESC = "Read summary statistics for case-control traits from the file specified";

	private final static String OPT_KEY = "key";
	private final static String OPT_KEY_DESC = "Self defined key workds: snp, beta, se, a1, a2, chr, bp, p";

	private final static String OPT_EXTRACT = "extract";
	private final static String OPT_EXTRACT_DESC = "Extract snps for joint mapping";
}
