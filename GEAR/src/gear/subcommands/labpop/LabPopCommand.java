package gear.subcommands.labpop;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class LabPopCommand extends Command
{

	@Override
	public String getName() 
	{
		return "labpop";
	}

	@Override
	public String getDescription()
	{
		return "Generating laborotary populations.";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options)
	{
		options.addOption(OptionBuilder.withDescription(OPT_SIZE_DESC).withLongOpt(OPT_SIZE_LONG).hasArg().create(OPT_SIZE));
		options.addOption(OptionBuilder.withDescription(OPT_NUM_MARKERS_DESC).withLongOpt(OPT_NUM_MARKERS_LONG).hasArg().create(OPT_NUM_MARKERS));
		options.addOption(OptionBuilder.withDescription(OPT_NULL_MARKER_LONG_DESC).withLongOpt(OPT_NULL_MARKER_LONG).hasArg().create());

		options.addOption(OptionBuilder.withDescription(OPT_SEED_DESC).withLongOpt(OPT_SEED_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_MAKE_BED_DESC).withLongOpt(OPT_MAKE_BED_LONG).create(OPT_MAKE_BED));

		options.addOption(OptionBuilder.withDescription(OPT_REC_DESC).withLongOpt(OPT_REC_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_RAND_REC_LONG_DESC).withLongOpt(OPT_RAND_REC_LONG).create());
		options.addOption(OptionBuilder.withDescription(OPT_REC_FILE_LONG_DESC).withLongOpt(OPT_REC_FILE_LONG).hasArg().create());

		options.addOption(OptionBuilder.withDescription(OPT_HSQ_DESC).hasArg().create(OPT_HSQ));
		options.addOption(OptionBuilder.withDescription(OPT_HSQ_DOM_DESC).hasArg().withLongOpt(OPT_HSQ_DOM).create());
		options.addOption(OptionBuilder.withDescription(OPT_HSQ_B_DESC).hasArg().withLongOpt(OPT_HSQ_B).create());

		options.addOption(OptionBuilder.withDescription(OPT_EFFECT_FILE_LONG_DESC).withLongOpt(OPT_EFFECT_FILE_LONG).hasArg().create());

		options.addOption(OptionBuilder.withDescription(OPT_DOM_EFFECT_FILE_LONG_DESC).withLongOpt(OPT_DOM_EFFECT_FILE_LONG).hasArg().create());

		options.addOption(OptionBuilder.withDescription(OPT_BC_DESC).create(OPT_BC));
		options.addOption(OptionBuilder.withDescription(OPT_F2_DESC).create(OPT_F2));
		options.addOption(OptionBuilder.withDescription(OPT_DH_DESC).create(OPT_DH));
		options.addOption(OptionBuilder.withDescription(OPT_RIL_DESC).create(OPT_RIL));
		options.addOption(OptionBuilder.withDescription(OPT_IF2_DESC).create(OPT_IF2));
		
		options.addOption(OptionBuilder.withDescription(OPT_REP_DESC).hasArg().create(OPT_REP));
		options.addOption(OptionBuilder.withDescription(OPT_EXCLUDE_DESC).hasArg().create(OPT_EXCLUDE));

		options.addOption(OptionBuilder.withDescription(OPT_ATGC_DESC).withLongOpt(OPT_ATGC_LONG).create());
		options.addOption(OptionBuilder.withDescription(OPT_1234_DESC).withLongOpt(OPT_1234_LONG).create());
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException
	{
		LabPopCommandArguments lpArgs = new LabPopCommandArguments();
		lpArgs.setSeed(parseLongOptionValue(cmdLine, OPT_SEED_LONG, OPT_SEED_DEFAULT));

		lpArgs.setSampleSize(parseIntOptionValue(cmdLine, OPT_SIZE_LONG, "100"));
		lpArgs.setNumberOfMarkers(parseIntOptionValue(cmdLine, OPT_NUM_MARKERS_LONG, "100"));
		lpArgs.setRec(parseDoubleOptionValueInRange(cmdLine, OPT_REC_LONG, "0.5", 0.0, 0.5));

		if (cmdLine.hasOption(OPT_NULL_MARKER_LONG))
		{
			lpArgs.setNullMarkerNum(cmdLine.getOptionValue(OPT_NULL_MARKER_LONG));
		}

//rec
		if (cmdLine.hasOption(OPT_REC_LONG))
		{
			lpArgs.setRec(parseDoubleOptionValueInRange(cmdLine, OPT_REC_LONG, "0.5", 0.0, 0.5));
		}
		if (cmdLine.hasOption(OPT_RAND_REC_LONG))
		{
			lpArgs.setRecRand();
		}
		if (cmdLine.hasOption(OPT_REC_FILE_LONG))
		{
			lpArgs.setRecFile(cmdLine.getOptionValue(OPT_REC_FILE_LONG));
		}
//ve
		if (cmdLine.hasOption(OPT_HSQ_B))
		{
			lpArgs.setHSQB(parseDoubleOptionValueInRange(cmdLine, OPT_HSQ_B, "0.5", 0.0, 0.99));
		}
		
//effect
		if (cmdLine.hasOption(OPT_HSQ))
		{
			lpArgs.setHsq(parseDoubleOptionValueInRange(cmdLine, OPT_HSQ, "0.5", 0.0, 0.99));
		}
		if (cmdLine.hasOption(OPT_EFFECT_FILE_LONG))
		{
			lpArgs.setPolyEffectFile(cmdLine.getOptionValue(OPT_EFFECT_FILE_LONG));
		}

//dom eff
		if (cmdLine.hasOption(OPT_DOM_EFFECT_FILE_LONG))
		{
			lpArgs.setPolyDomEffectFile(cmdLine.getOptionValue(OPT_DOM_EFFECT_FILE_LONG));
		}
		
		if (cmdLine.hasOption(OPT_HSQ_DOM))
		{
			lpArgs.setHsqDom(cmdLine.getOptionValue(OPT_HSQ_DOM));
		}
		
//pop
		if (cmdLine.hasOption(OPT_BC))
		{
			lpArgs.setBC();
		}

		if (cmdLine.hasOption(OPT_F2))
		{
			lpArgs.setF2();
		}

		if (cmdLine.hasOption(OPT_DH))
		{
			lpArgs.setDH();
		}

		if (cmdLine.hasOption(OPT_RIL))
		{
			lpArgs.setRIL();
		}
		
		if (cmdLine.hasOption(OPT_IF2))
		{
			lpArgs.setIF2();
		}

//bed		
		if (cmdLine.hasOption(OPT_MAKE_BED))
		{
			lpArgs.setMakeBed();
		}

		if (cmdLine.hasOption(OPT_REP))
		{
			lpArgs.setReplication(cmdLine.getOptionValue(OPT_REP));
		}

		if (cmdLine.hasOption(OPT_1234_LONG))
		{
			lpArgs.set1234mode();
		}
		if (cmdLine.hasOption(OPT_ATGC_LONG))
		{
			lpArgs.setATGCmode();
		}

		if (cmdLine.hasOption(OPT_EXCLUDE))
		{
			lpArgs.setExclude(cmdLine.getOptionValue(OPT_EXCLUDE));
		}
		return lpArgs;
	}

	@Override
	protected CommandImpl createCommandImpl() 
	{
		return new LabPopCommandImpl();
	}

	private static final char OPT_SIZE = 'n';
	private static final String OPT_SIZE_LONG = "sample-size";
	private static final String OPT_SIZE_DESC = "Specify the sample size.";

	private static final char OPT_NUM_MARKERS = 'm';
	private static final String OPT_NUM_MARKERS_LONG = "marker";
	private static final String OPT_NUM_MARKERS_DESC = "Specify the number of markers.";

	private static final String OPT_NULL_MARKER_LONG = "null-marker";
	private static final String OPT_NULL_MARKER_LONG_DESC = "Number of null markers.";

	private static final String OPT_REC_LONG = "rec";
	private static final String OPT_REC_DESC = "Specify the recombination fraction.";

	private static final String OPT_RAND_REC_LONG = "rand-rec";
	private static final String OPT_RAND_REC_LONG_DESC = "Use uniform distribution for recombination, between 0.01~0.5.";

	private static final String OPT_REC_FILE_LONG = "rec-file";
	private static final String OPT_REC_FILE_LONG_DESC = "Read recombination from the file specified. The first locus will converted to 0.5 regardless of value being specified.";

	
	//add-effect
	private static final String OPT_HSQ = "hsq";
	private static final String OPT_HSQ_DESC = "Heritability for polygenic model, 0.5 by default.";

	private static final String OPT_EFFECT_FILE_LONG = "effect-file";
	private static final String OPT_EFFECT_FILE_LONG_DESC = "Read effect from the file specified.";

	//dom-effect
	private static final String OPT_HSQ_DOM = "hsq-dom";
	private static final String OPT_HSQ_DOM_DESC = "heritability for dominance variance, 0 by default";

	private static final String OPT_DOM_EFFECT_FILE_LONG = "dom-effect-file";
	private static final String OPT_DOM_EFFECT_FILE_LONG_DESC = "Read dominance effect from the file specified.";

	private static final char OPT_MAKE_BED = 'b';
	private static final String OPT_MAKE_BED_LONG = "make-bed";
	private static final String OPT_MAKE_BED_DESC = "Make .bed, .bim and .fam files.";

	private static final String OPT_BC = "bc";
	private static final String OPT_BC_DESC = "Generate a backcross population.";

	private static final String OPT_F2 = "f2";
	private static final String OPT_F2_DESC = "Generate a F2 population.";

	private static final String OPT_DH = "dh";
	private static final String OPT_DH_DESC = "Generate a double haploid population.";

	private static final String OPT_RIL = "ril";
	private static final String OPT_RIL_DESC = "Generate a recombination inbred line population.";
	
	private static final String OPT_IF2 = "if2";
	private static final String OPT_IF2_DESC = "Generate an immortalized F2.";

	private static final String OPT_HSQ_B = "hsq-b";
	private static final String OPT_HSQ_B_DESC = "Broad-sense heritability";

	private static final String OPT_REP = "rep";
	private static final String OPT_REP_DESC = "Replication.";

	private static final String OPT_ATGC_LONG = "atgc-mode";
	private static final String OPT_ATGC_DESC = "Using atgc coding.";

	private static final String OPT_1234_LONG = "1234-mode";
	private static final String OPT_1234_DESC = "Using 1234 coding.";

	private static final String OPT_EXCLUDE = "exclude";
	private static final String OPT_EXCLUDE_DESC = "Exclude loci";
}
