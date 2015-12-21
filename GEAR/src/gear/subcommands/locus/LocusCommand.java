package gear.subcommands.locus;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.subcommands.eigengwas.EigenGWASArguments;

public class LocusCommand extends Command
{

	@Override
	public String getName()
	{
		return "locus";
	}

	@Override
	public String getDescription()
	{
		return "Statistics for loci";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options)
	{
		options.addOption(OptionBuilder.withDescription(OPT_BFILE_DESC).withLongOpt(OPT_BFILE_LONG).hasArg().isRequired().create());
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException
	{
		LocusArguments locusArgs = new LocusArguments();
		parseFileArguments(locusArgs, cmdLine);

		return locusArgs;
	}

	private void parseFileArguments(LocusArguments locusArgs, CommandLine cmdLine) throws CommandArgumentException
	{
		String bfile = cmdLine.getOptionValue("bfile");

		if ((bfile == null) )
		{
			throw new CommandArgumentException("No genotypes are provided.");
		}

		locusArgs.setBFile(bfile);
	}

	@Override
	protected CommandImpl createCommandImpl()
	{
		return new LocusImpl();
	}

}
