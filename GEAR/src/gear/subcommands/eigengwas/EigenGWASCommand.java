package gear.subcommands.eigengwas;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

public class EigenGWASCommand extends Command
{
  private static final String OPT_PHE_LONG = "pheno";
  private static final String OPT_PHE_DESC = "Specify the phenotype file (individual eigenvector)";
  private static final String OPT_MPHE_LONG = "mpheno";
  private static final String OPT_MPHE_DESC = "Specify the phenotype index";
  private static final String OPT_CHR = "chr";
  private static final String OPT_CHR_DESC = "Specify the chromosomes for analysis";

  public EigenGWASCommand()
  {
    addAlias("egwas");
  }

  public String getName()
  {
    return "eigengwas";
  }

  public String getDescription()
  {
    return "Eigen GWAS";
  }

  public void prepareOptions(Options options)
  {
    OptionBuilder.withDescription("Specify PLINK format .bed, .bim and .fam files"); OptionBuilder.withLongOpt("bfile"); OptionBuilder.hasArg(); OptionBuilder.isRequired(); options.addOption(OptionBuilder.create());
    OptionBuilder.withDescription("Specify the phenotype file (individual eigenvector)"); OptionBuilder.withLongOpt("pheno"); OptionBuilder.hasArg(); OptionBuilder.isRequired(); options.addOption(OptionBuilder.create());
    OptionBuilder.withDescription("Specify the phenotype index"); OptionBuilder.withLongOpt("mpheno"); OptionBuilder.hasArg(); options.addOption(OptionBuilder.create());
    OptionBuilder.withDescription("Specify the chromosomes for analysis"); OptionBuilder.hasArg(); options.addOption(OptionBuilder.create("chr"));
  }

  public CommandArguments parse(CommandLine cmdLine)
    throws CommandArgumentException
  {
    EigenGWASArguments eigenArgs = new EigenGWASArguments();
    parseFileArguments(eigenArgs, cmdLine);
    eigenArgs.setPhentypeIndex(parseIntOptionValue(cmdLine, "mpheno", "1"));
    eigenArgs.setPhenotypeFile(cmdLine.getOptionValue("pheno"));
    return eigenArgs;
  }

  private void parseFileArguments(EigenGWASArguments eigenArgs, CommandLine cmdLine) throws CommandArgumentException
  {
    String bfile = cmdLine.getOptionValue("bfile");
    String file = cmdLine.getOptionValue("file");

    if ((bfile == null) && (file == null))
    {
      throw new CommandArgumentException("No genotypes are provided. Either --bfile or --file must be set.");
    }

    if ((bfile != null) && (file != null))
    {
      throw new CommandArgumentException("--bfile and --file cannot be set together.");
    }

    eigenArgs.setBFile(bfile);
    eigenArgs.setFile(file);

    if (cmdLine.hasOption("chr"))
    {
      eigenArgs.setChr(cmdLine.getOptionValue("chr"));
    }
  }

  protected CommandImpl createCommandImpl()
  {
    return new EigenGWASImpl();
  }
}