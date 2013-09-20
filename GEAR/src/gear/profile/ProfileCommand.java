package gear.profile;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.Command;
import gear.CommandArgumentException;
import gear.CommandArguments;
import gear.CommandImpl;

public final class ProfileCommand extends Command
{

	@Override
	public String getName()
	{
		return "prof";
	}

	@Override
	public String getDescription()
	{
		return "Calculate the risk profile scores";
	}

	@SuppressWarnings("static-access")
	@Override
	protected void prepareOptions(Options options)
	{
		options.addOption(OptionBuilder.withDescription(OPT_SCORE_DESC).withLongOpt(OPT_SCORE_LONG).isRequired().hasArg().create(OPT_SCORE));
		options.addOption(OptionBuilder.withDescription(OPT_QSCORE_DESC).withLongOpt(OPT_QSCORE_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_QRANGE_DESC).withLongOpt(OPT_QRANGE_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_FILE_DESC).withLongOpt(OPT_FILE_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_BFILE_DESC).withLongOpt(OPT_BFILE_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_MACH_DOSAGE_DESC).withLongOpt(OPT_MACH_DOSAGE_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_MACH_INFO_DESC).withLongOpt(OPT_MACH_INFO_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_MACH_DOSAGE_BATCH_DESC).withLongOpt(OPT_MACH_DOSAGE_BATCH_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_MACH_INFO_BATCH_DESC).withLongOpt(OPT_MACH_INFO_BATCH_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_MODEL_DESC).withLongOpt(OPT_MODEL_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_LOGIT_DESC).withLongOpt(OPT_LOGIT_LONG).hasArg(false).create());
		options.addOption(OptionBuilder.withDescription(OPT_AUTO_FLIP_OFF_DESC).withLongOpt(OPT_AUTO_FLIP_OFF_LONG).hasArg(false).create());
		options.addOption(OptionBuilder.withDescription(OPT_NO_WEIGHT_DESC).withLongOpt(OPT_NO_WEIGHT_LONG).hasArg(false).create());
		options.addOption(OptionBuilder.withDescription(OPT_KEEP_ATGC_DESC).withLongOpt(OPT_KEEP_ATGC_LONG).hasArg(false).create());
		options.addOption(OptionBuilder.withDescription(OPT_OUT_DESC).withLongOpt(OPT_OUT_LONG).hasArg().create(OPT_OUT));
	}

	@Override
	protected CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException
	{
		ProfileCommandArguments profCmdArgs = new ProfileCommandArguments();
		profCmdArgs.setScoreFile(cmdLine.getOptionValue(OPT_SCORE_LONG));
		parseQScoreQRangeArgs(profCmdArgs, cmdLine);
		parseDataFileArgs(profCmdArgs, cmdLine);
		parseCoeffModelArgs(profCmdArgs, cmdLine);
		profCmdArgs.setIsLogit(cmdLine.hasOption(OPT_LOGIT_LONG));
		profCmdArgs.setIsAutoFlip(!cmdLine.hasOption(OPT_AUTO_FLIP_OFF_LONG));
		profCmdArgs.setIsWeighted(!cmdLine.hasOption(OPT_NO_WEIGHT_LONG));
		profCmdArgs.setIsKeepATGC(cmdLine.hasOption(OPT_KEEP_ATGC_LONG));
		profCmdArgs.setOutRoot(cmdLine.getOptionValue(OPT_OUT, OPT_OUT_DEFAULT));
		return profCmdArgs;
	}
	
	private void parseQScoreQRangeArgs(ProfileCommandArguments profCmdArgs, CommandLine cmdLine) throws CommandArgumentException
	{
		String qScoreFile = cmdLine.getOptionValue(OPT_QSCORE_LONG);
		String qRangeFile = cmdLine.getOptionValue(OPT_QRANGE_LONG);
		
		if (qScoreFile == null && qRangeFile != null || qScoreFile != null && qRangeFile == null)
		{
			throw new CommandArgumentException("--" + OPT_QSCORE_LONG + " and --" + OPT_QRANGE_LONG + " must be set together.");
		}
		
		profCmdArgs.setQScoreFile(qScoreFile);
		profCmdArgs.setQRangeFile(qRangeFile);
	}
	
	private void parseDataFileArgs(ProfileCommandArguments profCmdArgs, CommandLine cmdLine) throws CommandArgumentException
	{
		String file = cmdLine.getOptionValue(OPT_FILE_LONG);
		String bfile = cmdLine.getOptionValue(OPT_BFILE_LONG);
		String machDosageFile = cmdLine.getOptionValue(OPT_MACH_DOSAGE_LONG);
		String machInfoFile = cmdLine.getOptionValue(OPT_MACH_INFO_LONG);
		String machDosageBatch = cmdLine.getOptionValue(OPT_MACH_DOSAGE_BATCH_LONG);
		String machInfoBatch = cmdLine.getOptionValue(OPT_MACH_INFO_BATCH_LONG);
		
		if (file != null)
		{
			if (bfile != null)
			{
				throw new CommandArgumentException("--" + OPT_FILE_LONG + " cannot be set with --" + OPT_BFILE_LONG + ".");
			}
			else if (machDosageFile != null)
			{
				throw new CommandArgumentException("--" + OPT_FILE_LONG + " cannot be set with --" + OPT_MACH_DOSAGE_LONG + ".");
			}
			else if (machInfoFile != null)
			{
				throw new CommandArgumentException("--" + OPT_FILE_LONG + " cannot be set with --" + OPT_MACH_INFO_LONG + ".");
			}
			else if (machDosageBatch != null)
			{
				throw new CommandArgumentException("--" + OPT_FILE_LONG + " cannot be set with --" + OPT_MACH_DOSAGE_BATCH_LONG + ".");
			}
			else if (machInfoBatch != null)
			{
				throw new CommandArgumentException("--" + OPT_FILE_LONG + " cannot be set with --" + OPT_MACH_INFO_BATCH_LONG + ".");
			}
			profCmdArgs.setFile(file);
		}
		else if (bfile != null)
		{
			if (machDosageFile != null)
			{
				throw new CommandArgumentException("--" + OPT_BFILE_LONG + " cannot be set with --" + OPT_MACH_DOSAGE_LONG + ".");
			}
			else if (machInfoFile != null)
			{
				throw new CommandArgumentException("--" + OPT_BFILE_LONG + " cannot be set with --" + OPT_MACH_INFO_LONG + ".");
			}
			else if (machDosageBatch != null)
			{
				throw new CommandArgumentException("--" + OPT_BFILE_LONG + " cannot be set with --" + OPT_MACH_DOSAGE_BATCH_LONG + ".");
			}
			else if (machInfoBatch != null)
			{
				throw new CommandArgumentException("--" + OPT_BFILE_LONG + " cannot be set with --" + OPT_MACH_INFO_BATCH_LONG + ".");
			}
			profCmdArgs.setBFile(bfile);
		}
		else if (machDosageFile != null)
		{
			if (machInfoFile == null)
			{
				throw new CommandArgumentException("--" + OPT_MACH_DOSAGE_LONG + " must be set with --" + OPT_MACH_INFO_LONG + ".");
			}
			else if (machDosageBatch != null || machInfoBatch != null)
			{
				throw new CommandArgumentException("--" + OPT_MACH_DOSAGE_LONG + " cannot be set with --" + OPT_MACH_DOSAGE_BATCH_LONG + " or --" + OPT_MACH_INFO_BATCH_LONG + ".");
			}
			profCmdArgs.setMachDosageFile(machDosageFile);
			profCmdArgs.setMachInfoFile(machInfoFile);
		}
		else if (machInfoFile != null)
		{
			throw new CommandArgumentException("--" + OPT_MACH_INFO_LONG + " must be set with --" + OPT_MACH_DOSAGE_LONG + ".");
		}
		else if (machDosageBatch != null)
		{
			if (machInfoBatch == null)
			{
				throw new CommandArgumentException("--" + OPT_MACH_DOSAGE_BATCH_LONG + " must be set with --" + OPT_MACH_INFO_BATCH_LONG + ".");
			}
			profCmdArgs.setMachDosageBatch(machDosageBatch);
			profCmdArgs.setMachInfoBatch(machInfoBatch);
		}
		else if (machInfoBatch != null)
		{
			throw new CommandArgumentException("--" + OPT_MACH_INFO_BATCH_LONG + " must be set with --" + OPT_MACH_DOSAGE_BATCH_LONG + ".");
		}
		else
		{	
			throw new CommandArgumentException("No input data files.");
		}
	}
	
	private void parseCoeffModelArgs(ProfileCommandArguments profCmdArgs, CommandLine cmdLine) throws CommandArgumentException
	{
		String sCoeffModel = cmdLine.getOptionValue(OPT_MODEL_LONG);
		String resultFile = cmdLine.getOptionValue(OPT_OUT, OPT_OUT_DEFAULT);
		
		if (sCoeffModel == null)
		{
			profCmdArgs.setCoeffModel(new AdditiveCoeffModel());
			profCmdArgs.setIsSameAsPlink(true);
		}
		else
		{
			resultFile += "." + sCoeffModel;
			if (sCoeffModel.equals(OPT_MODEL_ADDITIVE))
			{
				profCmdArgs.setCoeffModel(new AdditiveCoeffModel());
			}
			else if (sCoeffModel.equals(OPT_MODEL_DOMINANCE))
			{
				profCmdArgs.setCoeffModel(new DominanceCoeffModel());
			}
			else if (sCoeffModel.equals(OPT_MODEL_RECESSIVE))
			{
				profCmdArgs.setCoeffModel(new RecessiveCoeffModel());
			}
			else
			{
				String msg = "";
				msg += "'" + sCoeffModel + "' is an invalid coefficient model. ";
				msg += "Valid models are '" + OPT_MODEL_ADDITIVE + "'(additive), ";
				msg += "'" + OPT_MODEL_DOMINANCE + "'(dominance) and ";
				msg += "'" + OPT_MODEL_RECESSIVE + "'(recessive).";
				throw new CommandArgumentException(msg);
			}
		}
		
		resultFile += ".profile";
		profCmdArgs.setResultFile(resultFile);
	}
	
	@Override
	protected CommandImpl createCommandImpl()
	{
		return new ProfileCommandImpl();
	}
	
	private static final char OPT_SCORE = 's';
	private static final String OPT_SCORE_LONG = "score";
	private static final String OPT_SCORE_DESC = "Specify score file";
	
	private static final String OPT_QSCORE_LONG = "qscore";
	private static final String OPT_QSCORE_DESC = "Specify q-score file";
	
	private static final String OPT_QRANGE_LONG = "qrange";
	private static final String OPT_QRANGE_DESC = "Specify q-score range file";
	
	private static final String OPT_MACH_DOSAGE_LONG = "mach-dosage";
	private static final String OPT_MACH_DOSAGE_DESC = "Specify MaCH format .mldose file";
	
	private static final String OPT_MACH_INFO_LONG = "mach-info";
	private static final String OPT_MACH_INFO_DESC = "Specify MaCH format .mlinfo file";
	
	private static final String OPT_MACH_DOSAGE_BATCH_LONG = "mach-dosage-batch";
	private static final String OPT_MACH_DOSAGE_BATCH_DESC = "Specify the text file in which each line records the name of a .mldose file";
	
	private static final String OPT_MACH_INFO_BATCH_LONG = "mach-info-batch";
	private static final String OPT_MACH_INFO_BATCH_DESC = "Specify the text file in which each line records the name of a .mlinfo file";
	
	private static final String OPT_MODEL_LONG = "model";
	private static final String OPT_MODEL_ADDITIVE = "add";
	private static final String OPT_MODEL_DOMINANCE = "dom";
	private static final String OPT_MODEL_RECESSIVE = "rec";
	private static final String OPT_MODEL_DESC = "Specify the compute Model, valid values are '" + OPT_MODEL_ADDITIVE + "'(additive), " +
	                                             "'" + OPT_MODEL_DOMINANCE + "'(dominance) and '" + OPT_MODEL_RECESSIVE + "'(recessive). " +
	                                             "If this option is not set, then PLINK model is used";
	
	private static final String OPT_LOGIT_LONG = "logit";
	private static final String OPT_LOGIT_DESC = "Take logit transform on the scores";
	
	private static final String OPT_AUTO_FLIP_OFF_LONG = "auto-flip-off";
	private static final String OPT_AUTO_FLIP_OFF_DESC = "Don't flip the alleles even if neither allele of a locus matches the score allele";
	
	private static final String OPT_NO_WEIGHT_LONG = "no-weight";
	private static final String OPT_NO_WEIGHT_DESC = "Don't divide the score sum by the number of loci";
	
	private static final String OPT_KEEP_ATGC_LONG = "keep-atgc";
	private static final String OPT_KEEP_ATGC_DESC = "Keep A/T and G/C loci";

}
