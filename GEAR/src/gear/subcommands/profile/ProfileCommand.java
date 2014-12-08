package gear.subcommands.profile;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.OptionGroup;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public final class ProfileCommand extends Command
{
	@Override
	public String getName()
	{
		return "profile";
	}

	@Override
	public String getDescription()
	{
		return "Calculate the risk profile scores";
	}
	
	public void setIsCalledByEnigma(boolean isCalledByEnigma)
	{
		this.isCalledByEnigma = isCalledByEnigma;
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options)
	{
		if (!isCalledByEnigma)
		{
			OptionGroup scoreFileGroup = new OptionGroup();
			scoreFileGroup.setRequired(true);
			scoreFileGroup.addOption(OptionBuilder.withDescription(OPT_SCORE_DESC).withLongOpt(OPT_SCORE_LONG).hasArg().create(OPT_SCORE));
			scoreFileGroup.addOption(OptionBuilder.withDescription(OPT_SCORE_GZ_DESC).withLongOpt(OPT_SCORE_GZ_LONG).hasArg().create());
			options.addOptionGroup(scoreFileGroup);
			options.addOption(OptionBuilder.withDescription(OPT_NO_SCORE_HEADER_DESC).withLongOpt(OPT_NO_SCORE_HEADER_LONG).create());
			options.addOption(OptionBuilder.withDescription(OPT_QSCORE_DESC).withLongOpt(OPT_QSCORE_LONG).hasArg().create());
			options.addOption(OptionBuilder.withDescription(OPT_QRANGE_DESC).withLongOpt(OPT_QRANGE_LONG).hasArg().create());
		}
		options.addOption(OptionBuilder.withDescription(OPT_FILE_DESC).withLongOpt(OPT_FILE_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_BFILE_DESC).withLongOpt(OPT_BFILE_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_MACH_DOSAGE_DESC).withLongOpt(OPT_MACH_DOSAGE_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_MACH_INFO_DESC).withLongOpt(OPT_MACH_INFO_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_MACH_DOSAGE_BATCH_DESC).withLongOpt(OPT_MACH_DOSAGE_BATCH_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_MACH_INFO_BATCH_DESC).withLongOpt(OPT_MACH_INFO_BATCH_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_MODEL_DESC).withLongOpt(OPT_MODEL_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_MODEL_FILE_DESC).withLongOpt(OPT_MODEL_FILE_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_LOGIT_DESC).withLongOpt(OPT_LOGIT_LONG).hasArg(false).create());
		options.addOption(OptionBuilder.withDescription(OPT_AUTO_FLIP_OFF_DESC).withLongOpt(OPT_AUTO_FLIP_OFF_LONG).hasArg(false).create());
		options.addOption(OptionBuilder.withDescription(OPT_NO_WEIGHT_DESC).withLongOpt(OPT_NO_WEIGHT_LONG).hasArg(false).create());
		options.addOption(OptionBuilder.withDescription(OPT_KEEP_ATGC_DESC).withLongOpt(OPT_KEEP_ATGC_LONG).hasArg(false).create());
		options.addOption(OptionBuilder.withDescription(OPT_EXTRACT_DESC).withLongOpt(OPT_EXTRACT_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_REMOVE_DESC).withLongOpt(OPT_REMOVE_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_SCALE_DESC).withLongOpt(OPT_SCALE_LONG).hasOptionalArg().create());
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException
	{
		ProfileCommandArguments profCmdArgs = new ProfileCommandArguments();
		profCmdArgs.setScoreFile(cmdLine.getOptionValue(OPT_SCORE_LONG));

		profCmdArgs.setScoreFileGZ(cmdLine.getOptionValue(OPT_SCORE_GZ_LONG));
		profCmdArgs.setHasScoreHeader(!cmdLine.hasOption(OPT_NO_SCORE_HEADER_LONG));
		parseQScoreQRangeArgs(profCmdArgs, cmdLine);
		parseDataFileArgs(profCmdArgs, cmdLine);
		parseCoeffModelArgs(profCmdArgs, cmdLine);
		profCmdArgs.setIsLogit(cmdLine.hasOption(OPT_LOGIT_LONG));
		profCmdArgs.setIsAutoFlip(!cmdLine.hasOption(OPT_AUTO_FLIP_OFF_LONG));
		profCmdArgs.setIsWeighted(!cmdLine.hasOption(OPT_NO_WEIGHT_LONG));
		profCmdArgs.setIsKeepATGC(cmdLine.hasOption(OPT_KEEP_ATGC_LONG));
		if (cmdLine.hasOption(OPT_EXTRACT_LONG))
		{
			profCmdArgs.setIsExtract(cmdLine.getOptionValue(OPT_EXTRACT_LONG));
		}

		if (cmdLine.hasOption(OPT_REMOVE_LONG))
		{
			profCmdArgs.setIsRemove(cmdLine.getOptionValue(OPT_REMOVE_LONG));
		}

		if (cmdLine.hasOption(OPT_SCALE_LONG))
		{
			profCmdArgs.setScale(parseStringOptionValue(cmdLine, OPT_SCALE_LONG, OPT_SCALE_DEFAULT));		
		}
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
		String coeffModel = cmdLine.getOptionValue(OPT_MODEL_LONG);
		String resultFile = cmdLine.getOptionValue(OPT_OUT, OPT_OUT_DEFAULT);
		
		profCmdArgs.setIsSameAsPlink(false);
		
		if (coeffModel == null)
		{
			String modelFile = cmdLine.getOptionValue(OPT_MODEL_FILE_LONG);
			if (modelFile == null)
			{
				profCmdArgs.setCoeffModelType(CoeffModelType.ADDITIVE);
				profCmdArgs.setIsSameAsPlink(true);
			}
			else
			{
				profCmdArgs.setCoeffModelType(CoeffModelType.FILE);
				profCmdArgs.setCoeffModelFile(modelFile);
			}
		}
		else
		{
			resultFile += "." + coeffModel;
			if (coeffModel.equals(OPT_MODEL_ADDITIVE))
			{
				profCmdArgs.setCoeffModelType(CoeffModelType.ADDITIVE);
			}
			else if (coeffModel.equals(OPT_MODEL_DOMINANCE))
			{
				profCmdArgs.setCoeffModelType(CoeffModelType.DOMINANCE);
			}
			else if (coeffModel.equals(OPT_MODEL_RECESSIVE))
			{
				profCmdArgs.setCoeffModelType(CoeffModelType.RECESSIVE);
			}
			else
			{
				String msg = "";
				msg += "'" + coeffModel + "' is an invalid coefficient model. ";
				msg += "Valid models are '" + OPT_MODEL_ADDITIVE + "'(additive), ";
				msg += "'" + OPT_MODEL_DOMINANCE + "'(dominance) and ";
				msg += "'" + OPT_MODEL_RECESSIVE + "'(recessive).";
				throw new CommandArgumentException(msg);
			}
		}
		
		profCmdArgs.setResultFile(resultFile);
	}
	
	@Override
	protected CommandImpl createCommandImpl()
	{
		return new ProfileCommandImpl();
	}

	private boolean isCalledByEnigma;
	
	private static final char OPT_SCORE = 's';
	private static final String OPT_SCORE_LONG = "score";
	private static final String OPT_SCORE_DESC = "Specify score file";
	
	private static final String OPT_SCORE_GZ_LONG = "score-gz";
	private static final String OPT_SCORE_GZ_DESC = "Specify score file in gzip format";
	
	private static final String OPT_NO_SCORE_HEADER_LONG = "no-score-header";
	private static final String OPT_NO_SCORE_HEADER_DESC = "Indicate that the score file has no header (title) line";
	
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
	
	private static final String OPT_MODEL_FILE_LONG = "model-file";
	private static final String OPT_MODEL_FILE_DESC = "Specify the model file";
	
	private static final String OPT_MODEL_LONG = "model";
	private static final String OPT_MODEL_ADDITIVE = "add";
	private static final String OPT_MODEL_DOMINANCE = "dom";
	private static final String OPT_MODEL_RECESSIVE = "rec";
	private static final String OPT_MODEL_DESC = "Specify the compute Model, valid values are '" + OPT_MODEL_ADDITIVE + "'(additive), " +
	                                             "'" + OPT_MODEL_DOMINANCE + "'(dominance) and '" + OPT_MODEL_RECESSIVE + "'(recessive). " +
	                                             "If this option and --" + OPT_MODEL_FILE_LONG + " is not set, then PLINK model is used";
	
	private static final String OPT_LOGIT_LONG = "logit";
	private static final String OPT_LOGIT_DESC = "Take logit transform on the scores";
	
	private static final String OPT_AUTO_FLIP_OFF_LONG = "auto-flip-off";
	private static final String OPT_AUTO_FLIP_OFF_DESC = "Don't flip the alleles even if neither allele of a locus matches the score allele";
	
	private static final String OPT_NO_WEIGHT_LONG = "no-weight";
	private static final String OPT_NO_WEIGHT_DESC = "Don't divide the score sum by the number of loci";
	
	private static final String OPT_KEEP_ATGC_LONG = "keep-atgc";
	private static final String OPT_KEEP_ATGC_DESC = "Keep A/T and G/C loci";
	
	private static final String OPT_EXTRACT_LONG = "extract-score";
	private static final String OPT_EXTRACT_DESC = "Extract score snps";

	private static final String OPT_REMOVE_LONG = "remove-score";
	private static final String OPT_REMOVE_DESC = "Remove score snps";
	
	private static final String OPT_SCALE_LONG = "scale";
	private static final String OPT_SCALE_DEFAULT = null;
	private static final String OPT_SCALE_DESC = "Standardise genotypes";
}
