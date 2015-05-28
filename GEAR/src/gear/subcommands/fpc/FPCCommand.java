package gear.subcommands.fpc;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class FPCCommand extends Command
{

	@Override
	public String getName()
	{
		return "fpc";
	}

	@Override
	public String getDescription()
	{
		return "Inference of Fst-derived principal components.";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options)
	{
//		options.addOption(OptionBuilder.withDescription(OPT_META_BATCH_DESC).withLongOpt(OPT_META_BATCH_LONG).hasArg().create(OPT_META_BATCH));
		options.addOption(OptionBuilder.withDescription(OPT_FST_DESC).hasArg().create(OPT_FST));
		options.addOption(OptionBuilder.withDescription(OPT_REF_DESC).hasArgs(3).create(OPT_REF));
		options.addOption(OptionBuilder.withDescription(OPT_COORDINATE_LONG_DESC).withLongOpt(OPT_COORDINATE_LONG).hasArg().create(OPT_COORDINATE));
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine)
			throws CommandArgumentException
	{
		FPCCommandArguments FPCArg = new FPCCommandArguments();
		FPCArg.setFstFile(cmdLine.getOptionValue(OPT_FST));
//		FPCArg.setMetaBatch(cmdLine.getOptionValue(OPT_META_BATCH));

		if (cmdLine.hasOption(OPT_REF))
		{
			FPCArg.setReference(cmdLine.getOptionValues(OPT_REF));
		}

		if (cmdLine.hasOption(OPT_COORDINATE))
		{
			FPCArg.setCoordinates(cmdLine.getOptionValue(OPT_COORDINATE));
		}
		return FPCArg;
	}

	@Override
	protected CommandImpl createCommandImpl()
	{
		return new FPCCommandImpl();
	}

//	private final static String OPT_META_BATCH = "mb";
//	private final static String OPT_META_BATCH_LONG = "meta-batch";
//	private final static String OPT_META_BATCH_DESC = "The summary statistic batch";
//
	private final static String OPT_FST = "fst";
	private final static String OPT_FST_DESC = "Specify the fst file.";

	private final static String OPT_REF = "ref";
	private final static String OPT_REF_DESC = "Specify the index for reference populations.";

	private final static String OPT_COORDINATE = "c";
	private final static String OPT_COORDINATE_LONG = "cordinates";
	private final static String OPT_COORDINATE_LONG_DESC = "Specify the coordinates.";

}
