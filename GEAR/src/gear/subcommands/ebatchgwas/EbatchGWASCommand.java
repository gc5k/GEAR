package gear.subcommands.ebatchgwas;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class EbatchGWASCommand extends Command
{

	@Override
	public String getName()
	{
		return "eigengwas";
	}

	@Override
	public String getDescription()
	{
		return "all in one for eigengwas";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options)
	{
		options.addOption(OptionBuilder.withDescription(OPT_BFILE_DESC).withLongOpt(OPT_BFILE_LONG).hasArg().isRequired().create());
		options.addOption(OptionBuilder.withDescription(OPT_EV_DESC).withLongOpt(OPT_EV_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_DOM_DESC).create(OPT_DOM));

	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException
	{
		EbatchGWASArguments EArgs = new EbatchGWASArguments();
		parseFileArguments(EArgs, cmdLine);
		EArgs.setEV(cmdLine.getOptionValue(OPT_EV_LONG));
		EArgs.setDom(cmdLine.hasOption(OPT_DOM));
		return EArgs;
	}

	private void parseFileArguments(EbatchGWASArguments EArgs, CommandLine cmdLine) throws CommandArgumentException
	{
		String bfile = cmdLine.getOptionValue("bfile");

		if ((bfile == null) )
		{
			throw new CommandArgumentException("No genotypes are provided. Either --bfile or --file must be set.");
		}

		EArgs.setBFile(bfile);
	}

	@Override
	protected CommandImpl createCommandImpl() 
	{
		return new EbatchGWASImpl();
	}

	private final static String OPT_EV_LONG = "ev";
	private final static String OPT_EV_DESC = "Specify the eigenvector numbers";
	private final static String OPT_DOM = "dom";
	private final static String OPT_DOM_DESC = "Specify the eigenvector numbers";
}
