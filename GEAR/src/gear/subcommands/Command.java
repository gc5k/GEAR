package gear.subcommands;

import gear.AboutInfo;
import gear.util.FileUtil;
import gear.util.Logger;

import java.util.Collections;
import java.util.Iterator;
import java.util.Set;
import java.util.TreeSet;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;

public abstract class Command implements Comparable<Command> {
	public abstract String getName();

	public boolean hasAlias(String alias) {
		return aliases.contains(alias);
	}

	protected void addAlias(String alias) {
		aliases.add(alias);
	}

	public Set<String> getAliases() {
		return Collections.unmodifiableSet(aliases);
	}

	@Override
	public boolean equals(Object obj) {
		String s = null;

		if (obj instanceof String) {
			s = (String) obj;
		} else if (obj instanceof Command) {
			s = ((Command) obj).getName();
		}

		return s != null && (getName().equals(s) || hasAlias(s));
	}

	public int compareTo(Command otherCmd) {
		if (equals(otherCmd)) {
			return 0;
		}
		return getName().compareTo(otherCmd.getName());
	}

	private Set<String> aliases = new TreeSet<String>();

	public abstract String getDescription();

	public String getLongDescription() {
		return getDescription();
	}

	public String getFullDescription() {
		return "";
	}

	public abstract void prepareOptions(Options options);

	@SuppressWarnings("static-access")
	public Options getOptions() {
		Options options = new Options();
		options.addOption(
				OptionBuilder.withDescription(OPT_OUT_DESC).withLongOpt(OPT_OUT_LONG).hasArg().create());
		prepareOptions(options);
		return options;
	}

	public abstract CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException;

	protected abstract CommandImpl createCommandImpl();

	protected void parseMAFArguments(CommandArguments cmdArgs, CommandLine cmdLine) throws CommandArgumentException {
		if (cmdLine.hasOption(OPT_MAF_LONG)) {
			cmdArgs.setMAF(cmdLine.getOptionValue(OPT_MAF_LONG));
		}
	}

	protected void parseMAXMAFArguments(CommandArguments cmdArgs, CommandLine cmdLine) throws CommandArgumentException {
		if (cmdLine.hasOption(OPT_MAX_MAF_LONG)) {
			cmdArgs.setMaxMAF(cmdLine.getOptionValue(OPT_MAX_MAF_LONG));
		}
	}

	protected void parseMAFRangeArguments(CommandArguments cmdArgs, CommandLine cmdLine) throws CommandArgumentException {
		if (cmdLine.hasOption(OPT_MAF_RANGE_LONG)) {
			String[] mafR = cmdLine.getOptionValues(OPT_MAF_RANGE_LONG);
			for (int i = 0; i < mafR.length; i++) {
				if(!mafR[i].contains("-")) {
					Logger.printUserLog("The value of " + mafR[i] + " is incorrect for option --" + OPT_MAF_RANGE_LONG + ". GEAR quit.");
					System.exit(1);
				}
			}
			cmdArgs.setMAFRange(cmdLine.getOptionValues(OPT_MAF_RANGE_LONG));
		}
	}

	protected void parseGENOArguments(CommandArguments cmdArgs, CommandLine cmdLine) throws CommandArgumentException {
		if (cmdLine.hasOption(OPT_GENO_LONG)) {
			cmdArgs.setMaxMAF(cmdLine.getOptionValue(OPT_GENO_LONG));
		}
	}

	protected void parseZeroVarArguments(CommandArguments cmdArgs, CommandLine cmdLine) throws CommandArgumentException {
		if (cmdLine.hasOption(OPT_ZERO_VAR_LONG)) {
			cmdArgs.setZeroVar();
		}
	}

	protected void parseFileArguments(CommandArguments cmdArgs, CommandLine cmdLine) throws CommandArgumentException {
		String bfile = cmdLine.getOptionValue(OPT_BFILE_LONG);
		String file = cmdLine.getOptionValue(OPT_FILE_LONG);

		if ((bfile == null) && (file == null)) {
			throw new CommandArgumentException("No genotypes are provided. Either --bfile or --file must be set.");
		}

		if ((bfile != null) && (file != null)) {
			throw new CommandArgumentException("--bfile and --file cannot be set together.");
		}

		if (bfile != null) {
			FileUtil.exists(new String(bfile + ".bed"));
			FileUtil.exists(new String(bfile + ".bim"));
			FileUtil.exists(new String(bfile + ".fam"));
			cmdArgs.setBFile(bfile);
		}
		if (file != null) {
			FileUtil.exists(new String(bfile + ".ped"));
			FileUtil.exists(new String(bfile + ".map"));
			cmdArgs.setFile(file);
		}
	}

	protected void parsePhenoFileArguments(CommandArguments cmdArgs, CommandLine cmdLine) 
			throws CommandArgumentException {
		if (cmdLine.hasOption(OPT_PHE)) {
			FileUtil.exists(cmdLine.getOptionValue(OPT_PHE));
			cmdArgs.setPhenotypeFile(cmdLine.getOptionValue(OPT_PHE));
		}
	}

	protected void parsePhenoIndexArguments(CommandArguments cmdArgs, CommandLine cmdLine) 
			throws CommandArgumentException {
		if (cmdLine.hasOption(OPT_MPHE)) {
			cmdArgs.setPhenotypeIndex(cmdLine.getOptionValues(OPT_MPHE));
		}
	}

	protected void parseSampleFilterArguments(CommandArguments cmdArgs, CommandLine cmdLine)
			throws CommandArgumentException {
		String keepFile = cmdLine.getOptionValue(OPT_KEEP_LONG);
		String removeFile = cmdLine.getOptionValue(OPT_REMOVE_LONG);

		if (cmdLine.hasOption(OPT_KEEP_LONG) && cmdLine.hasOption(OPT_REMOVE_LONG)) {
			throw new CommandArgumentException("--" + OPT_KEEP_LONG +" and --" + OPT_KEEP_LONG +" cannot be set together.");
		}
		if (keepFile != null) {
			FileUtil.exists(keepFile);
			cmdArgs.setKeepFile(keepFile);
		}
		if (removeFile != null) {
			FileUtil.exists(removeFile);
			cmdArgs.setRemoveFile(removeFile);
		}
	}

	protected void parseSNPFilterFileArguments(CommandArguments cmdArgs, CommandLine cmdLine)
		throws CommandArgumentException {

		if (cmdLine.hasOption(OPT_EXTRACT_LONG) && cmdLine.hasOption(OPT_EXCLUDE_LONG)) {
			throw new CommandArgumentException(
				"--extract and --exclude cannot be set together.");
		}

		if (cmdLine.hasOption(OPT_EXTRACT_LONG)) {
			FileUtil.exists(cmdLine.getOptionValue(OPT_EXTRACT_LONG));
			cmdArgs.setExtractFile(cmdLine.getOptionValue(OPT_EXTRACT_LONG));
		} else if (cmdLine.hasOption(OPT_EXCLUDE_LONG)) {
			FileUtil.exists(cmdLine.getOptionValue(OPT_EXCLUDE_LONG));
			cmdArgs.setExcludeFile(cmdLine.getOptionValue(OPT_EXCLUDE_LONG));
		}
	}

	protected void parseSNPFilterChromosomeArguments(CommandArguments cmdArgs, CommandLine cmdLine)
			throws CommandArgumentException {
		if (cmdLine.hasOption(OPT_CHR_LONG)
				&& cmdLine.hasOption(OPT_NOT_CHR_LONG)) {
			throw new CommandArgumentException(
					"--chr and --not-chr cannot be set together.");
		}

		if (cmdLine.hasOption(OPT_CHR_LONG)) {
			cmdArgs.setChr(cmdLine.getOptionValues(OPT_CHR_LONG));
		} else if (cmdLine.hasOption(OPT_NOT_CHR_LONG)) {
			cmdArgs.setNotChr(cmdLine.getOptionValues(OPT_NOT_CHR_LONG));
		}
	}

	protected void printOptionsInEffect(CommandLine cmdLine, String subcmd) {
		Logger.printUserLog("Subcommand: " + subcmd);
		Logger.printUserLog("Options in effect: ");

		@SuppressWarnings("rawtypes")
		Iterator optIter = cmdLine.iterator();
		while (optIter.hasNext()) {
			Option opt = (Option) optIter.next();
			String line = opt.hasLongOpt() ? "\t--" + opt.getLongOpt() : "\t--" + opt.getOpt();
			String[] argValues = opt.getValues();
			if (argValues != null) {
				for (String value : argValues) {
					line += " " + value;
				}
			}
			Logger.printUserLog(line);
		}

		Logger.printUserLog("");
	}

	public void execute(String[] args, String subCmdName) {
		Options options = getOptions();
		CommandLineParser cmdLineParser = new PosixParser();

		try {
			CommandLine cmdLine = cmdLineParser.parse(options, args, stopAtNonOption);
			CommandArguments cmdArgs = parse(cmdLine);
			cmdArgs.setOutRoot(cmdLine.getOptionValue(OPT_OUT_LONG, OPT_OUT_DEFAULT));

			Logger.setLogFiles(cmdArgs.getOutRoot());
			Logger.hasUserLogTag(false);
			Logger.printUserLog(AboutInfo.WELCOME_MESSAGE);
			Logger.hasUserLogTag(true);

			printOptionsInEffect(cmdLine, subCmdName);

			CommandImpl cmdImpl = createCommandImpl();
			cmdImpl.preExecute();
			cmdImpl.execute(cmdArgs);
			cmdImpl.postExecute();
		} catch (ParseException e) {
			Logger.printUserError(e.getMessage());
			System.exit(1);
		} catch (CommandArgumentException e) {
			Logger.printUserError(e.getMessage());
			System.exit(1);
		}
	}

	protected void setStopAtNonOption(boolean stopAtNonOption) {
		this.stopAtNonOption = stopAtNonOption;
	}

	protected int parseIntOptionValue(CommandLine cmdLine, String opt, String defaultOptVal)
			throws CommandArgumentException {
		int value;
		try {
			value = Integer.parseInt(cmdLine.getOptionValue(opt, defaultOptVal));
		} catch (NumberFormatException e) {
			String msg = "";
			msg += "The value of --" + opt + "is invalid: '";
			msg += cmdLine.getOptionValue(opt) + "' is not a valid integer.";
			throw new CommandArgumentException(msg);
		}
		return value;
	}

	protected long parseLongOptionValue(CommandLine cmdLine, String opt) throws CommandArgumentException {
		return parseLongOptionValue(cmdLine, opt, null);
	}

	protected long parseLongOptionValue(CommandLine cmdLine, String opt, String defaultOptVal)
			throws CommandArgumentException {
		long value;
		try {
			value = Long.parseLong(cmdLine.getOptionValue(opt, defaultOptVal));
		} catch (NumberFormatException e) {
			String msg = "";
			msg += "The value of --" + opt + "is invalid: '";
			msg += cmdLine.getOptionValue(opt) + "' is not a valid integer.";
			throw new CommandArgumentException(msg);
		}
		return value;
	}

	protected long parseLongOptionValueInRange(CommandLine cmdLine, String opt, String defaultOptVal, long min,
			long max) throws CommandArgumentException {
		long value = parseLongOptionValue(cmdLine, opt, defaultOptVal);
		if (value < min || value > max) {
			throw new CommandArgumentException(
					"--" + opt + " must be no smaller than " + min + " and no larger than " + max);
		}
		return value;
	}

	protected double parseDoubleOptionValue(CommandLine cmdLine, String opt) throws CommandArgumentException {
		return parseDoubleOptionValue(cmdLine, opt, null);
	}

	protected double parseDoubleOptionValue(CommandLine cmdLine, String opt, String defaultOptVal)
			throws CommandArgumentException {
		double value;
		try {
			value = Double.parseDouble(cmdLine.getOptionValue(opt, defaultOptVal));
		} catch (NumberFormatException e) {
			String msg = "";
			msg += "The value of --" + opt + "is invalid: '";
			msg += cmdLine.getOptionValue(opt) + "' is not a valid floating point number.";
			throw new CommandArgumentException(msg);
		}
		return value;
	}

	protected double parseDoubleOptionValueInRange(CommandLine cmdLine, String opt, String defaultOptVal, double min,
			double max) throws CommandArgumentException {
		double value = parseDoubleOptionValue(cmdLine, opt, defaultOptVal);
		if (value < min || value > max) {
			throw new CommandArgumentException(
					"--" + opt + " must be no smaller than " + min + " and no larger than " + max);
		}
		return value;
	}

	protected String parseStringOptionValue(CommandLine cmdLine, String opt, String defaultOptVal)
			throws CommandArgumentException {
		String value = cmdLine.getOptionValue(opt, defaultOptVal);
		return value;
	}

	private boolean stopAtNonOption;

	protected static final String OPT_MAF_LONG = "maf";
	protected static final String OPT_MAF_DESC = "Minimal allele frequency cutoff. By dafault 0.01.";

	protected static final String OPT_MAX_MAF_LONG = "max-maf";
	protected static final String OPT_MAX_MAF_DESC = "Upper bound for minimal allele frequency cutoff. By dafault 0.5.";

	protected static final String OPT_MAF_RANGE_LONG = "maf-range";
	protected static final String OPT_MAF_RANGE_DESC = "maf ranges.";

	protected static final String OPT_GENO_LONG = "genotyping missing rate";
	protected static final String OPT_GENO_DESC = "Genotyping missing rate cutoff. By default 0.1";

	protected static final String OPT_ZERO_VAR_LONG = "zero-var";
	protected static final String OPT_ZERO_VAR_DESC = "Zero variance loci.";

	protected static final String OPT_KEEP_LONG = "keep";
	protected static final String OPT_KEEP_DESC = "Specify the individuals for analysis";

	protected static final String OPT_REMOVE_LONG = "remove";
	protected static final String OPT_REMOVE_DESC = "Specify the individuals removed from analysis";

	protected static final String OPT_CHR_LONG = "chr";
	protected static final String OPT_CHR_DESC = "Specify the chromosomes for analysis";

	protected static final String OPT_NOT_CHR_LONG = "not-chr";
	protected static final String OPT_NOT_CHR_DESC = "Specify the chromosomes not for analysis";

	protected static final String OPT_EXTRACT_LONG = "extract";
	protected static final String OPT_EXTRACT_DESC = "Specify the SNPs for analysis";

	protected static final String OPT_EXCLUDE_LONG = "exclude";
	protected static final String OPT_EXCLUDE_DESC = "Specify the SNPs excluded from analysis";

	protected static final String OPT_FILE_LONG = "file";
	protected static final String OPT_FILE_DESC = "Specify PLINK format .ped and .map files";

	protected static final String OPT_BFILE_LONG = "bfile";
	protected static final String OPT_BFILE_DESC = "Specify PLINK format .bed, .bim and .fam files";

	protected static final String OPT_PHE = "pheno";
	protected static final String OPT_PHE_DESC = "Specify the phenotype file (individual eigenvector)";

	protected static final String OPT_MPHE = "mpheno";
	protected static final String OPT_MPHE_DESC = "Specify the phenotype index";

	protected static final char OPT_OUT = 'o';
	protected static final String OPT_OUT_LONG = "out";
	protected static final String OPT_OUT_DEFAULT = "gear";
	protected static final String OPT_OUT_DESC = "Specify output root filename, default to '" + OPT_OUT_DEFAULT + "'";

	protected static final String OPT_SEED_LONG = "seed";
	protected static final String OPT_SEED_DEFAULT = "2012";
	protected static final String OPT_SEED_DESC = "Specify the seed of random number generator, default to "
			+ OPT_SEED_DEFAULT;
}
