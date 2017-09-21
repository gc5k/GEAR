package gear.subcommands.simulationfam;

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
		options.addOption(OptionBuilder.withDescription(OPT_FREQ_FILE_LONG_DESC).withLongOpt(OPT_FREQ_FILE_LONG).hasArg().create());

		options.addOption(OptionBuilder.withDescription(OPT_REC_DESC).withLongOpt(OPT_REC_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_REC_SEX_DESC).withLongOpt(OPT_REC_SEX_LONG).hasArgs(2).create(OPT_REC_SEX));
		options.addOption(OptionBuilder.withDescription(OPT_REC_UNIF_DESC).withLongOpt(OPT_REC_UNIF_LONG).create());

//		options.addOption(OptionBuilder.withDescription(OPT_QTL_DESC).withLongOpt(OPT_QTL_LONG).hasArg().create(OPT_QTL));
		options.addOption(OptionBuilder.withDescription(OPT_HSQ_DESC).hasArg().create(OPT_HSQ));
		options.addOption(OptionBuilder.withDescription(OPT_HSQ_DOM_DESC).hasArg().withLongOpt(OPT_HSQ_DOM).create());
//		options.addOption(OptionBuilder.withDescription("Scale IBD").create(OPT_SCALE));
		options.addOption(OptionBuilder.withDescription(OPT_EFFECT_FILE_LONG_DESC).withLongOpt(OPT_EFFECT_FILE_LONG).hasArg().create());

		options.addOption(OptionBuilder.withDescription(OPT_DOM_EFFECT_DESC).hasArg().withLongOpt(OPT_DOM_EFFECT).create());
		options.addOption(OptionBuilder.withDescription(OPT_POLY_DOM_EFFECT_LONG_DESC).withLongOpt(OPT_POLY_DOM_EFFECT_LONG).create());
		options.addOption(OptionBuilder.withDescription(OPT_POLY_DOM_EFFECT_SORT_LONG_DESC).withLongOpt(OPT_POLY_DOM_EFFECT_SORT_LONG).create());
		options.addOption(OptionBuilder.withDescription(OPT_DOM_EFFECT_FILE_LONG_DESC).withLongOpt(OPT_DOM_EFFECT_FILE_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_REP_DESC).hasArg().create(OPT_REP));
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException
	{
		SimuFamilyCommandArguments famArgs = new SimuFamilyCommandArguments();
		famArgs.setNumberOfFamilies(parseIntOptionValue(cmdLine, OPT_NUM_FAMS_LONG, "100"));
		famArgs.setNumberOfMarkers(parseIntOptionValue(cmdLine, OPT_NUM_MARKERS_LONG, "100"));
		if(cmdLine.hasOption(OPT_NULL_MARKER_LONG))
		{
			famArgs.setNullMarkerNum(cmdLine.getOptionValue(OPT_NULL_MARKER_LONG));
		}

		famArgs.setSeed(parseLongOptionValue(cmdLine, OPT_SEED_LONG, OPT_SEED_DEFAULT));

//maf
		famArgs.setMAF(parseDoubleOptionValueInRange(cmdLine, OPT_MAF, "0.5", 0.01, 0.500001));
		if (cmdLine.hasOption(OPT_UNIF_MAF_LONG))
		{
			famArgs.setUnifMAF();
		}

		if (cmdLine.hasOption(OPT_FREQ_FILE_LONG))
		{
			famArgs.setFreqFile(cmdLine.getOptionValue(OPT_FREQ_FILE_LONG));
		}

//ld
		famArgs.setLD(parseDoubleOptionValueInRange(cmdLine, OPT_LD_LONG, "0", -1, 1));
		if (cmdLine.hasOption(OPT_RAND_LD_LONG))
		{
			famArgs.setRandLD();
		}

//rec
		famArgs.setRec(parseDoubleOptionValueInRange(cmdLine, OPT_REC_LONG, "0.5", 0.01, 0.500001));
		if(cmdLine.hasOption(OPT_REC_SEX_LONG))
		{
			famArgs.setRecSex(cmdLine.getOptionValues(OPT_REC_SEX_LONG));
		}
		if(cmdLine.hasOption(OPT_REC_UNIF_LONG))
		{
			famArgs.setRecRandFlag();
		}

//h2
//		if(cmdLine.hasOption(OPT_QTL))
//		{
//			cmdArgs.setQTLFile(cmdLine.getOptionValue(OPT_QTL));
//		}
		if(cmdLine.hasOption(OPT_HSQ))
		{
			famArgs.setHsq(cmdLine.getOptionValue(OPT_HSQ));
		}

		if(cmdLine.hasOption(OPT_EFFECT_FILE_LONG))
		{
			famArgs.setPolyEffectFile(cmdLine.getOptionValue(OPT_EFFECT_FILE_LONG));
		}

		if(cmdLine.hasOption(OPT_MAKE_BED))
		{
			famArgs.setMakeBed();
		}

//		if(cmdLine.hasOption(OPT_SCALE))
//		{
//			famArgs.setScale();
//		}
		
		//dom eff
		famArgs.setPolyDomEffect();

		if (cmdLine.hasOption(OPT_DOM_EFFECT))
		{
			famArgs.setPlainDomEffect(parseDoubleOptionValue(cmdLine, OPT_DOM_EFFECT, "0.5"));
		}

		if (cmdLine.hasOption(OPT_POLY_DOM_EFFECT_SORT_LONG))
		{
			famArgs.setPolyDomEffectSort();
		}

		if (cmdLine.hasOption(OPT_DOM_EFFECT_FILE_LONG))
		{
			famArgs.setPolyDomEffectFile(cmdLine.getOptionValue(OPT_DOM_EFFECT_FILE_LONG));
		}

		if (cmdLine.hasOption(OPT_HSQ_DOM))
		{
			famArgs.setHsqDom(cmdLine.getOptionValue(OPT_HSQ_DOM));
		}

		if (cmdLine.hasOption(OPT_REP))
		{
			famArgs.setRep(cmdLine.getOptionValue(OPT_REP));
		}
		return famArgs;
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

	private static final String OPT_FREQ_FILE_LONG = "freq-file";
	private static final String OPT_FREQ_FILE_LONG_DESC = "Read frequency from the file specified.";

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

	private static final String OPT_EFFECT_FILE_LONG = "effect-file";
	private static final String OPT_EFFECT_FILE_LONG_DESC = "Read effect from the file specified.";

	//dom-effect
	private static final String OPT_HSQ_DOM = "hsq-dom";
	private static final String OPT_HSQ_DOM_DESC = "heritability for dominance variance";
	private static final String OPT_DOM_EFFECT = "dom-effect";
	private static final String OPT_DOM_EFFECT_DESC = "Equal dominant effects, 1 by default.";

	private static final String OPT_POLY_DOM_EFFECT_LONG = "poly-dom-effect";
	private static final String OPT_POLY_DOM_EFFECT_LONG_DESC = "Polygenic dominance (normal distribution) model.";

	private static final String OPT_POLY_DOM_EFFECT_SORT_LONG = "poly-effect-sort";
	private static final String OPT_POLY_DOM_EFFECT_SORT_LONG_DESC = "Sorted polygenic (normal distribution) model.";

	private static final String OPT_DOM_EFFECT_FILE_LONG = "dom-effect-file";
	private static final String OPT_DOM_EFFECT_FILE_LONG_DESC = "Read dominance effect from the file specified.";
	
	private static final String OPT_REP = "rep";
	private static final String OPT_REP_DESC = "replication.";

}
