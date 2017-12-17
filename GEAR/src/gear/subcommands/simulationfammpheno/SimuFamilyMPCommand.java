package gear.subcommands.simulationfammpheno;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public final class SimuFamilyMPCommand extends Command
{
	@Override
	public String getName()
	{
		return "simufammp";
	}

	@Override
	public String getDescription()
	{
		return "Simulation of discordant nuclear families";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options)
	{
		options.addOption(OptionBuilder.withDescription(OPT_NUM_FAMS_DESC).withLongOpt(OPT_NUM_FAMS_LONG).hasArg().create(OPT_NUM_FAMS));
		options.addOption(OptionBuilder.withDescription(OPT_NUM_MARKERS_DESC).withLongOpt(OPT_NUM_MARKERS_LONG).hasArg().create(OPT_NUM_MARKERS));
		options.addOption(OptionBuilder.withDescription(OPT_NULL_MARKER_LONG_DESC).withLongOpt(OPT_NULL_MARKER_LONG).hasArg().create());

		options.addOption(OptionBuilder.withDescription(OPT_SEED_DESC).withLongOpt(OPT_SEED_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_MAKE_BED_DESC).withLongOpt(OPT_MAKE_BED_LONG).create(OPT_MAKE_BED));

		options.addOption(OptionBuilder.withDescription(OPT_LD_DESC).withLongOpt(OPT_LD_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_RAND_LD_DESC).withLongOpt(OPT_RAND_LD_LONG).create());

		options.addOption(OptionBuilder.withDescription(OPT_MAF_DESC).hasArg().create(OPT_MAF));
		options.addOption(OptionBuilder.withDescription(OPT_UNIF_MAF_LONG_DESC).withLongOpt(OPT_UNIF_MAF_LONG).create());
		options.addOption(OptionBuilder.withDescription(OPT_REC_DESC).withLongOpt(OPT_REC_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_REC_SEX_DESC).withLongOpt(OPT_REC_SEX_LONG).hasArgs(2).create(OPT_REC_SEX));
		options.addOption(OptionBuilder.withDescription(OPT_REC_UNIF_DESC).withLongOpt(OPT_REC_UNIF_LONG).create());

//		options.addOption(OptionBuilder.withDescription(OPT_QTL_DESC).withLongOpt(OPT_QTL_LONG).hasArg().create(OPT_QTL));
		options.addOption(OptionBuilder.withDescription(OPT_HSQ_DESC).hasArgs().create(OPT_HSQ));
		options.addOption(OptionBuilder.withDescription(OPT_COR_MAT_LONG_DESC).withLongOpt(OPT_COR_MAT_LONG).hasArg().create(OPT_CM));
		options.addOption(OptionBuilder.withDescription(OPT_COR_MAT_E_LONG_DESC).withLongOpt(OPT_COR_MAT_E_LONG).hasArg().create(OPT_CME));
		
		options.addOption(OptionBuilder.withDescription(OPT_EFFECT_FILE_LONG_DESC).withLongOpt(OPT_EFFECT_FILE_LONG).hasArg().create());
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException
	{
		SimuFamilyMPCommandArguments famMPArgs = new SimuFamilyMPCommandArguments();
		famMPArgs.setNumberOfFamilies(parseIntOptionValue(cmdLine, OPT_NUM_FAMS_LONG, "100"));
		famMPArgs.setNumberOfMarkers(parseIntOptionValue(cmdLine, OPT_NUM_MARKERS_LONG, "100"));
		if(cmdLine.hasOption(OPT_NULL_MARKER_LONG))
		{
			famMPArgs.setNullMarkerNum(cmdLine.getOptionValue(OPT_NULL_MARKER_LONG));
		}

		famMPArgs.setSeed(parseLongOptionValue(cmdLine, OPT_SEED_LONG, OPT_SEED_DEFAULT));

//maf
		famMPArgs.setMAF(parseDoubleOptionValueInRange(cmdLine, OPT_MAF, "0.5", 0.01, 0.500001));
		if (cmdLine.hasOption(OPT_UNIF_MAF_LONG))
		{
			famMPArgs.setUnifMAF();
		}

//ld
		famMPArgs.setLD(parseDoubleOptionValueInRange(cmdLine, OPT_LD_LONG, "0", -1, 1));			
		if (cmdLine.hasOption(OPT_RAND_LD_LONG))
		{
			famMPArgs.setRandLD();
		}

//rec
		famMPArgs.setRec(parseDoubleOptionValueInRange(cmdLine, OPT_REC_LONG, "0.5", 0.01, 0.500001));
		if(cmdLine.hasOption(OPT_REC_SEX_LONG))
		{
			famMPArgs.setRecSex(cmdLine.getOptionValues(OPT_REC_SEX_LONG));
		}
		if(cmdLine.hasOption(OPT_REC_UNIF_LONG))
		{
			famMPArgs.setRecRandFlag();
		}

//h2
//		if(cmdLine.hasOption(OPT_QTL))
//		{
//			cmdArgs.setQTLFile(cmdLine.getOptionValue(OPT_QTL));
//		}
		if(cmdLine.hasOption(OPT_HSQ))
		{
			famMPArgs.setHsq(cmdLine.getOptionValues(OPT_HSQ));
		}

		if(cmdLine.hasOption(OPT_EFFECT_FILE_LONG))
		{
			famMPArgs.setPolyEffectFile(cmdLine.getOptionValue(OPT_EFFECT_FILE_LONG));
		}

		if (cmdLine.hasOption(OPT_CM))
		{
			famMPArgs.setCM(cmdLine.getOptionValue(OPT_CM));
		}

		if (cmdLine.hasOption(OPT_CME))
		{
			famMPArgs.setCME(cmdLine.getOptionValue(OPT_CME));
		}

		if(cmdLine.hasOption(OPT_MAKE_BED))
		{
			famMPArgs.setMakeBed();
		}
		return famMPArgs;
	}

	@Override
	protected CommandImpl createCommandImpl()
	{
		return new SimuFamilyMPCommandImpl();
	}
	
	private static final char OPT_NUM_FAMS = 'f';
	private static final String OPT_NUM_FAMS_LONG = "num-fam";
	private static final String OPT_NUM_FAMS_DESC = "Specify the number of families";
	
	private static final char OPT_NUM_MARKERS = 'm';
	private static final String OPT_NUM_MARKERS_LONG = "num-marker";
	private static final String OPT_NUM_MARKERS_DESC = "Specify the number of markers";
	
	private static final String OPT_NULL_MARKER_LONG = "null-marker";
	private static final String OPT_NULL_MARKER_LONG_DESC = "Number of null markers.";

	private static final char OPT_MAKE_BED = 'b';
	private static final String OPT_MAKE_BED_LONG = "make-bed";
	private static final String OPT_MAKE_BED_DESC = "Make .bed, .bim and .fam files";

	private static final String OPT_LD_LONG = "ld";
	private static final String OPT_LD_DESC = "Specify the ld (Lewontin's DPrime)";

	private static final String OPT_RAND_LD_LONG = "rand-ld";
	private static final String OPT_RAND_LD_DESC = "Generate the ld (Lewontin's DPrime) from uniform distribtuion between -1 and 1";

	private static final String OPT_MAF = "freq";
	private static final String OPT_MAF_DESC = "Specify the allele frequency";

	private static final String OPT_UNIF_MAF_LONG = "unif-freq";
	private static final String OPT_UNIF_MAF_LONG_DESC = "Use uniform distribution for MAF, between 0.01~0.5.";

	private static final String OPT_REC_LONG = "rec";
	private static final String OPT_REC_DESC = "Specify the recombination fraction";

	private static final String OPT_REC_SEX = "rs";
	private static final String OPT_REC_SEX_LONG = "rec-sex";
	private static final String OPT_REC_SEX_DESC = "Specify the sex-specific recombination fraction";

	private static final String OPT_REC_UNIF_LONG = "unif-rec";
	private static final String OPT_REC_UNIF_DESC = "Use uniform distribution recombination fractions beween (0~0.5)";

//	private static final String OPT_QTL = "q";
//	private static final String OPT_QTL_LONG = "qtl";
//	private static final String OPT_QTL_DESC = "qtl parameters (locp, locm, effm, effp, h2)";
	
	private static final String OPT_HSQ = "hsq";
	private static final String OPT_HSQ_DESC = "Heritability for polygenic model, 0.5 by default.";
	
	private static final String OPT_COR_MAT_LONG = "cor-mat";
	private static final String OPT_CM = "cm";
	private static final String OPT_COR_MAT_LONG_DESC = "correlation matrix";

	private static final String OPT_COR_MAT_E_LONG = "cor-mat-e";
	private static final String OPT_CME = "cme";
	private static final String OPT_COR_MAT_E_LONG_DESC = "correlation matrix for environment";

	private static final String OPT_EFFECT_FILE_LONG = "effect-file";
	private static final String OPT_EFFECT_FILE_LONG_DESC = "Read effect from the file specified.";


}
