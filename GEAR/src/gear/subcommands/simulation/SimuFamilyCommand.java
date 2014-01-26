package gear.subcommands.simulation;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public final class SimuFamilyCommand extends Command
{

	@Override
	public String getName()
	{
		return "simufam";
	}

	@Override
	public String getDescription()
	{
		return "Simulation of discordant nuclear families";
	}

	@SuppressWarnings("static-access")
	@Override
	protected void prepareOptions(Options options)
	{
		options.addOption(OptionBuilder.withDescription(OPT_NUM_FAMS_DESC).withLongOpt(OPT_NUM_FAMS_LONG).hasArg().isRequired().create(OPT_NUM_FAMS));
		options.addOption(OptionBuilder.withDescription(OPT_NUM_MARKERS_DESC).withLongOpt(OPT_NUM_MARKERS_LONG).hasArg().isRequired().create(OPT_NUM_MARKERS));
		options.addOption(OptionBuilder.withDescription(OPT_SEED_DESC).withLongOpt(OPT_SEED_LONG).hasArg().create(OPT_SEED));
		options.addOption(OptionBuilder.withDescription(OPT_MAKE_BED_DESC).withLongOpt(OPT_MAKE_BED_LONG).create(OPT_MAKE_BED));

		options.addOption(OptionBuilder.withDescription(OPT_LD_DESC).withLongOpt(OPT_LD_LONG).hasArg().create(OPT_LD));
		options.addOption(OptionBuilder.withDescription(OPT_REC_DESC).withLongOpt(OPT_REC_LONG).hasArg().create(OPT_REC));
		options.addOption(OptionBuilder.withDescription(OPT_REC_SEX_DESC).withLongOpt(OPT_REC_SEX_LONG).hasArg().create(OPT_REC_SEX));
		options.addOption(OptionBuilder.withDescription(OPT_REC_RAND_DESC).withLongOpt(OPT_REC_RAND_LONG).create(OPT_REC_RAND));

		options.addOption(OptionBuilder.withDescription(OPT_QTL_DESC).withLongOpt(OPT_QTL_LONG).hasArg().create(OPT_QTL));
	}

	@Override
	protected CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException
	{
		SimuFamilyCommandArguments cmdArgs = new SimuFamilyCommandArguments();
		parseNumberOfFamilies(cmdArgs, cmdLine);
		parseNumberOfMarkers(cmdArgs, cmdLine);
		parseSeed(cmdArgs, cmdLine);
		parseLD(cmdArgs, cmdLine);
		parseRec(cmdArgs, cmdLine);
		parseRecSex(cmdArgs, cmdLine);
		parseQTL(cmdArgs, cmdLine);
		cmdArgs.setRecRandFlag(cmdLine.hasOption(OPT_REC_RAND));
		cmdArgs.setRecSexFlag(cmdLine.hasOption(OPT_REC_SEX));
		cmdArgs.setMakeBed(cmdLine.hasOption(OPT_MAKE_BED));
		return cmdArgs;
	}

	private void parseNumberOfFamilies(SimuFamilyCommandArguments cmdArgs, CommandLine cmdLine) throws CommandArgumentException
	{
		int numFams = 0;
		try
		{
			numFams = Integer.parseInt(cmdLine.getOptionValue(OPT_NUM_FAMS));
		}
		catch (NumberFormatException e)
		{
		}
		
		if (numFams <= 0)
		{
			String msg = "";
			msg += "The value of --" + OPT_NUM_FAMS_LONG + "is invalid: '";
			msg += cmdLine.getOptionValue(OPT_NUM_FAMS) + "' is not a valid positive integer.";
			throw new CommandArgumentException(msg);
		}
		
		cmdArgs.setNumberOfFamilies(numFams);
	}
	
	private void parseNumberOfMarkers(SimuFamilyCommandArguments cmdArgs, CommandLine cmdLine) throws CommandArgumentException
	{
		int numMarkers = 0;
		
		try
		{
			numMarkers = Integer.parseInt(cmdLine.getOptionValue(OPT_NUM_MARKERS));
		}
		catch (NumberFormatException e)
		{
		}
		
		if (numMarkers <= 0)
		{
			String msg = "";
			msg += "The value of --" + OPT_NUM_MARKERS_LONG + "is invalid: '";
			msg += cmdLine.getOptionValue(OPT_NUM_MARKERS) + "' is not a valid positive integer.";
			throw new CommandArgumentException(msg);
		}
		
		cmdArgs.setNumberOfMarkers(numMarkers);
	}
	
	private void parseSeed(SimuFamilyCommandArguments cmdArgs, CommandLine cmdLine) throws CommandArgumentException
	{
		if (cmdLine.hasOption(OPT_SEED))
		{
			try
			{
				cmdArgs.setSeed(Long.parseLong(cmdLine.getOptionValue(OPT_SEED)));
			}
			catch (NumberFormatException e)
			{
				String msg = "";
				msg += "The value of --" + OPT_SEED_LONG + "is invalid: '";
				msg += cmdLine.getOptionValue(OPT_SEED) + "' is not a valid integer.";
				throw new CommandArgumentException(msg);
			}
		}
	}

	private void parseRec(SimuFamilyCommandArguments cmdArgs, CommandLine cmdLine) throws CommandArgumentException
	{
		double r = 0.5;
		if(cmdLine.hasOption(OPT_REC))
		{
			try
			{
				r = Double.parseDouble(cmdLine.getOptionValue(OPT_REC));
			}
			catch (NumberFormatException e)
			{
			}

			if (r < 0 || r > 0.5)
			{
				String msg = "";
				msg += "The value of --" + OPT_REC_LONG + "is invalid: '";
				msg += cmdLine.getOptionValue(OPT_REC) + "' is not a valid number between 0 and 0.5.";
				throw new CommandArgumentException(msg);
			}			
		}
		cmdArgs.setRec(r);
	}

	private void parseRecSex(SimuFamilyCommandArguments cmdArgs, CommandLine cmdLine) throws CommandArgumentException
	{
		double[] rs = {0.5, 0.5};
		if(cmdLine.hasOption(OPT_REC_SEX))
		{
			try
			{
				String[] s = cmdLine.getOptionValue(OPT_REC_SEX).split(",");
				rs[0] = Double.parseDouble(s[0]);
				rs[1] = Double.parseDouble(s[1]);
			}
			catch (NumberFormatException e)
			{
			}

			if (rs[0] < 0 || rs[0] > 0.5)
			{
				String msg = "";
				msg += "The value of --" + OPT_REC_SEX_LONG + "is invalid: '";
				msg += cmdLine.getOptionValue(OPT_REC_SEX) + "' is not a valid number between 0 and 0.5.";
				throw new CommandArgumentException(msg);
			}
			
			if (rs[1] < 0 || rs[1] > 0.5)
			{
				String msg = "";
				msg += "The value of --" + OPT_REC_SEX_LONG + "is invalid: '";
				msg += cmdLine.getOptionValue(OPT_REC_SEX) + "' is not a valid number between 0 and 0.5.";
				throw new CommandArgumentException(msg);
			}
		}
		cmdArgs.setRecSex(rs);
	}

	private void parseQTL(SimuFamilyCommandArguments cmdArgs, CommandLine cmdLine) throws CommandArgumentException
	{
		String f = null;
		
		if (cmdLine.hasOption(OPT_QTL))
		{
			f = cmdLine.getOptionValue(OPT_QTL);
		}
		cmdArgs.setQTLFile(f);
	}

	private void parseLD(SimuFamilyCommandArguments cmdArgs, CommandLine cmdLine) throws CommandArgumentException
	{
		double l = 0;
		
		if (cmdLine.hasOption(OPT_LD))
		{
			try
			{
				l = Double.parseDouble(cmdLine.getOptionValue(OPT_LD));
			}
			catch (NumberFormatException e)
			{
			}

			if (l < 0 || l > 1)
			{
				String msg = "";
				msg += "The value of --" + OPT_LD + "is invalid: '";
				msg += cmdLine.getOptionValue(OPT_LD) + "' is not a valid between 0 and 1.";
				throw new CommandArgumentException(msg);
			}			
		}
		cmdArgs.setLD(l);
	}

	@Override
	protected CommandImpl createCommandImpl()
	{
		return new SimuFamilyCommandImpl();
	}
	
	private static final char OPT_NUM_FAMS = 'f';
	private static final String OPT_NUM_FAMS_LONG = "num-fam";
	private static final String OPT_NUM_FAMS_DESC = "Specify the number of families";
	
	private static final char OPT_NUM_MARKERS = 'm';
	private static final String OPT_NUM_MARKERS_LONG = "num-marker";
	private static final String OPT_NUM_MARKERS_DESC = "Specify the number of markers";
	
	private static final char OPT_SEED = 's';
	private static final String OPT_SEED_LONG = "seed";
	private static final String OPT_SEED_DESC = "Specify the seed of random number generator";
	
	private static final char OPT_MAKE_BED = 'b';
	private static final String OPT_MAKE_BED_LONG = "make-bed";
	private static final String OPT_MAKE_BED_DESC = "Make .bed, .bim and .fam files";

	private static final char OPT_LD = 'l';
	private static final String OPT_LD_LONG = "ld";
	private static final String OPT_LD_DESC = "Specify the ld (Lewontin's DPrime)";

	private static final char OPT_REC = 'r';
	private static final String OPT_REC_LONG = "rec";
	private static final String OPT_REC_DESC = "Specify the recombination fraction";

	private static final String OPT_REC_SEX = "rs";
	private static final String OPT_REC_SEX_LONG = "rec-sex";
	private static final String OPT_REC_SEX_DESC = "Specify the sex-specific recombination fraction";

	private static final String OPT_REC_RAND = "rr";
	private static final String OPT_REC_RAND_LONG = "rec-rand";
	private static final String OPT_REC_RAND_DESC = "Use uniformly distributed recombination fractions beween (0~0.5)";

	private static final String OPT_QTL = "q";
	private static final String OPT_QTL_LONG = "qtl";
	private static final String OPT_QTL_DESC = "qtl parameters (locp, locm, effm, effp, h2)";
}
