package gear.subcommands.ppcbatch;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.OptionGroup;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.subcommands.exsnp.ExSNPCommand;
import gear.subcommands.exsnp.ExSNPCommandArguments;
import gear.subcommands.profile.ProfileCommand;
import gear.subcommands.profile.ProfileCommandArguments;

public class PPCBatchCommand extends Command
{

	@Override
	public String getName()
	{
		return "probatch";
	}

	@Override
	public String getDescription()
	{
		return "Generate projected PC from the batch file.";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options)
	{
		OptionGroup baseFileGroup = new OptionGroup();
		baseFileGroup.setRequired(true);

		baseFileGroup.addOption(OptionBuilder.withDescription(OPT_SCORE_DESC).withLongOpt(OPT_SCORE_LONG).hasArg().create(OPT_SCORE));
		baseFileGroup.addOption(OptionBuilder.withDescription(OPT_SCORE_GZ_DESC).withLongOpt(OPT_SCORE_GZ_LONG).hasArg().create());
		options.addOptionGroup(baseFileGroup);

		options.addOption(OptionBuilder.withDescription(OPT_BATCH_DESC).hasArg().create(OPT_BATCH));
		options.addOption(OptionBuilder.withDescription(OPT_NO_SCORE_HEADER_DESC).withLongOpt(OPT_NO_SCORE_HEADER_LONG).create());
//		options.addOption(OptionBuilder.withDescription(OPT_QSCORE_DESC).withLongOpt(OPT_QSCORE_LONG).hasArg().create());
//		options.addOption(OptionBuilder.withDescription(OPT_QRANGE_DESC).withLongOpt(OPT_QRANGE_LONG).hasArg().create());

		options.addOption(OptionBuilder.withDescription(OPT_GREEDY_DESC).create(OPT_GREEDY));
		exsnpCommand.prepareOptions(options);
		profileCommand.setIsCalledByEnigma(true);
		profileCommand.prepareOptions(options);
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine)
			throws CommandArgumentException
	{
		PPCBatchCommandArguments pbCmdArgs = new PPCBatchCommandArguments();
		
		pbCmdArgs.setBatch(cmdLine.getOptionValue(OPT_BATCH));
		
		if (cmdLine.hasOption(OPT_NO_SCORE_HEADER_LONG))
		{
			pbCmdArgs.setHeader();
		}

		pbCmdArgs.setScoreFile(cmdLine.getOptionValue(OPT_SCORE_LONG));
		pbCmdArgs.setScoreFileGZ(cmdLine.getOptionValue(OPT_SCORE_GZ_LONG));
		
		if (cmdLine.hasOption(OPT_GREEDY))
		{
			pbCmdArgs.setGreedy();
		}
		
		pbCmdArgs.setExSNPCommandArguments((ExSNPCommandArguments) exsnpCommand.parse(cmdLine));
		pbCmdArgs.setProfileCommandArguments((ProfileCommandArguments) profileCommand.parse(cmdLine));
		return pbCmdArgs;
	}

	@Override
	protected CommandImpl createCommandImpl()
	{
		return new PPCBatchCommandImpl();
	}

	private static final String OPT_BATCH = "batch";
	private static final String OPT_BATCH_DESC = "Specify batch file";

	private static final char OPT_SCORE = 's';
	private static final String OPT_SCORE_LONG = "score";
	private static final String OPT_SCORE_DESC = "Specify score file";

	private static final String OPT_SCORE_GZ_LONG = "score-gz";
	private static final String OPT_SCORE_GZ_DESC = "Specify score file in gzip format";

	private static final String OPT_NO_SCORE_HEADER_LONG = "no-score-header";
	private static final String OPT_NO_SCORE_HEADER_DESC = "Indicate that the score file has no header (title) line";

//	private static final String OPT_QSCORE_LONG = "qscore";
//	private static final String OPT_QSCORE_DESC = "Specify q-score file";

//	private static final String OPT_QRANGE_LONG = "qrange";
//	private static final String OPT_QRANGE_DESC = "Specify q-score range file";
	
	private static final String OPT_GREEDY = "greedy";
	private static final String OPT_GREEDY_DESC = "Using all markers in each bed file.";

	private ExSNPCommand exsnpCommand = new ExSNPCommand();
	private ProfileCommand profileCommand = new ProfileCommand();
}
