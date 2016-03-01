package gear.subcommands.at;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class AtCommand extends Command
{

	public AtCommand()
	{
		addAlias("ld");
	}

	@Override
	public String getName() 
	{
		return "ld-score";
	}
	@Override
	public String getDescription() 
	{
		// TODO Auto-generated method stub
		return "Generate ld scores for chromosomes";
	}
	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options)
	{
		options.addOption(OptionBuilder.withDescription(OPT_BFILE_DESC).withLongOpt(OPT_BFILE_LONG).hasArg().isRequired().create());
		options.addOption(OptionBuilder.withDescription(OPT_DPRIME_DESC).create(OPT_DPRIME));
		options.addOption(OptionBuilder.withDescription(OPT_D_DESC).create(OPT_D));
		options.addOption(OptionBuilder.withDescription(OPT_RSQ_DESC).create(OPT_RSQ));
		options.addOption(OptionBuilder.withDescription(OPT_WINDOW_DESC).hasArg().create(OPT_WINDOW));

	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException {
		AtCommandArguments atArgs = new AtCommandArguments();
		parseFileArguments(atArgs, cmdLine);

		if(cmdLine.hasOption(OPT_WINDOW))
		{
			atArgs.setWindow(cmdLine.getOptionValue(OPT_WINDOW));
		}

		if(cmdLine.hasOption(OPT_DPRIME))
		{
			atArgs.setDPrime();
		}

		if(cmdLine.hasOption(OPT_D))
		{
			atArgs.setD();
		}

		if(cmdLine.hasOption(OPT_RSQ))
		{
			atArgs.setRsq();
		}

		return atArgs;
	}
	
	private void parseFileArguments(AtCommandArguments atArgs, CommandLine cmdLine) throws CommandArgumentException
	{
		String bfile = cmdLine.getOptionValue("bfile");

		if (bfile == null)
		{
			throw new CommandArgumentException("No genotypes are provided. Either --bfile or --file must be set.");
		}

		atArgs.setBFile(bfile);
	}
	
	@Override
	protected CommandImpl createCommandImpl()
	{
		return new AtCommandImpl();
	}

	private static final String OPT_WINDOW = "window";
	private static final String OPT_WINDOW_DESC = "SNP pairs having distance smaller than the values (mb) specified will have ld calculation";

	private static final String OPT_DPRIME = "dprime";
	private static final String OPT_DPRIME_DESC = "Lewontin's Dprime, Genetics, 1964, 49:49-67";

	private static final String OPT_D = "d";
	private static final String OPT_D_DESC = "D for LD measure, d=correlation*sqrt(p1q1p2q2).";

	private static final String OPT_RSQ = "rsq";
	private static final String OPT_RSQ_DESC = "Rsq for LD measure.";

}
