package gear.subcommands.fst;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class FstCommand extends Command {
    
    @Override
    public String getName() {
        return "fst";
    }
    
    @Override
    public String getDescription() {
        return "Estimating fst";
    }
    
    @SuppressWarnings("static-access")
    @Override
    public void prepareOptions(Options options) {
        options.addOption(OptionBuilder.withDescription(OPT_BFILE_DESC).withLongOpt(OPT_BFILE_LONG).hasArg().isRequired().create());
        options.addOption(OptionBuilder.withDescription(OPT_GROUP_DESC).withLongOpt(OPT_GROUP_LONG).hasArg().isRequired().create());
        
        options.addOption(OptionBuilder.withDescription(OPT_KEEP_DESC).withLongOpt(OPT_KEEP_LONG).hasArg().create());
        options.addOption(
                          OptionBuilder.withDescription(OPT_REMOVE_DESC).withLongOpt(OPT_REMOVE_LONG).hasArg().create());
        options.addOption(OptionBuilder.withDescription(OPT_KEEP_FAM_DESC).withLongOpt(OPT_KEEP_FAM_LONG).hasArg().create());
        options.addOption(
                          OptionBuilder.withDescription(OPT_REMOVE_FAM_DESC).withLongOpt(OPT_REMOVE_FAM_LONG).hasArg().create());
        
        options.addOption(
                          OptionBuilder.withDescription(OPT_EXTRACT_DESC).withLongOpt(OPT_EXTRACT_LONG).hasArg().create());
        options.addOption(
                          OptionBuilder.withDescription(OPT_EXCLUDE_DESC).withLongOpt(OPT_EXCLUDE_LONG).hasArg().create());
        
        options.addOption(OptionBuilder.withDescription(OPT_CHR_DESC).withLongOpt(OPT_CHR_LONG).hasArgs().create());
        options.addOption(
                          OptionBuilder.withDescription(OPT_NOT_CHR_DESC).withLongOpt(OPT_NOT_CHR_LONG).hasArgs().create());
        
        options.addOption(OptionBuilder.withDescription(OPT_MAF_DESC).withLongOpt(OPT_MAF_LONG).hasArg().create());
        options.addOption(
                          OptionBuilder.withDescription(OPT_MAX_MAF_DESC).withLongOpt(OPT_MAX_MAF_LONG).hasArg().create());
        options.addOption(OptionBuilder.withDescription(OPT_GENO_DESC).withLongOpt(OPT_GENO_LONG).hasArg().create());
        options.addOption(
                          OptionBuilder.withDescription(OPT_ZERO_VAR_DESC).withLongOpt(OPT_ZERO_VAR_LONG).create());
        options.addOption(
                          OptionBuilder.withDescription(OPT_MAF_RANGE_DESC).withLongOpt(OPT_MAF_RANGE_LONG).hasArgs().create());
        
    }
    
    @Override
    public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException {
        FSTCommandArguments fstArgs = new FSTCommandArguments();
        
        parseFileArguments((CommandArguments) fstArgs, cmdLine);
        parseSampleFilterArguments((CommandArguments) fstArgs, cmdLine);
        parseFamilyFilterArguments((CommandArguments) fstArgs, cmdLine);
        
        parseSNPFilterFileArguments((CommandArguments) fstArgs, cmdLine);
        parseSNPFilterChromosomeArguments((CommandArguments) fstArgs, cmdLine);
        
        parseMAFArguments((CommandArguments) fstArgs, cmdLine);
        parseMAXMAFArguments((CommandArguments) fstArgs, cmdLine);
        parseGENOArguments((CommandArguments) fstArgs, cmdLine);
        parseZeroVarArguments((CommandArguments) fstArgs, cmdLine);
        
        parseMAFRangeArguments((CommandArguments) fstArgs, cmdLine);
        
        fstArgs.setGroup(cmdLine.getOptionValue(OPT_GROUP_LONG));
        return fstArgs;
    }
    
    @Override
    protected CommandImpl createCommandImpl() {
        return new FSTCommandImpl();
    }
    
    private static String OPT_GROUP_LONG = "group";
    private static String OPT_GROUP_DESC = "group";
    
}
