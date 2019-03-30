package gear;

//import java.util.Calendar;
import java.util.Iterator;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

//import gear.he.h2trans.H2Transformer;
//import gear.pscontrol.NonTransmitted;
import gear.subcommands.Command;
//import gear.util.Logger;
//import gear.util.MonitorThread;

public enum Gear {
	INSTANCE;
	public static String subcmdName;

	private Gear() {
		// hello, git!
		addCommand(new gear.subcommands.bluppca.BlupPcaCommand());
		addCommand(new gear.subcommands.dnafingerprint.DFPCommand());
		addCommand(new gear.subcommands.he.assocpower.HEAssocPowerCommand());
		addCommand(new gear.subcommands.help.HelpCommand());
		addCommand(new gear.subcommands.hpc.HpcCommand());
		addCommand(new gear.subcommands.lambdaD.LambdaDCommand());
		addCommand(new gear.subcommands.metahet.MetaHetCommand());
		addCommand(new gear.subcommands.weightedmeta.WeightedMetaCommand());
		addCommand(new gear.subcommands.mlmmeta.MLMMetaCommand());
		addCommand(new gear.subcommands.glsmeta.GLSMetaCommand());
		addCommand(new gear.subcommands.metawatchdog.powercalculator.DogPowerCommand());
		addCommand(new gear.subcommands.metawatchdog.decrypt.MetaWatchdog2Command());
		addCommand(new gear.subcommands.metawatchdog.encrypt.EnigmaCommand());
		addCommand(new gear.subcommands.profile.ProfileCommand());
		addCommand(new gear.subcommands.exsnp.ExSNPCommand());
		addCommand(new gear.subcommands.exsnp2.ExSNP2Command());
		addCommand(new gear.subcommands.lmdr.LMDRCommand());
		addCommand(new gear.subcommands.ppcbatch.PPCBatchCommand());
		addCommand(new gear.subcommands.projectedpc.ProjectedPCCommand());
		addCommand(new gear.subcommands.fst.FstCommand());

		addCommand(new gear.subcommands.sfst.SFstCommand());
		addCommand(new gear.subcommands.fpc.FPCCommand());
		addCommand(new gear.subcommands.metapc.MetaPCCommand());

		addCommand(new gear.subcommands.grmstat.GRMStatCommand());
		addCommand(new gear.subcommands.wgrmA.WGRMACommand());

		addCommand(new gear.subcommands.qpca.QPCACommand());
		addCommand(new gear.subcommands.locus.LocusCommand());
		addCommand(new gear.subcommands.locusA.LocusACommand());

		addCommand(new gear.subcommands.at.AtCommand());
		addCommand(new gear.subcommands.ebatchgwas.EbatchGWASCommand());
		addCommand(new gear.subcommands.eigengwas.EigenGWASCommand());
		addCommand(new gear.subcommands.eigengwasdom.EigenGWASDomCommand());
		addCommand(new gear.subcommands.eigengwasepi.EigenGWASEpiCommand());

		addCommand(new gear.subcommands.labpop.LabPopCommand());
		addCommand(new gear.subcommands.oath.nss.NSSCommand());
		addCommand(new gear.subcommands.oath.oathx.OATHXCommand());
		addCommand(new gear.subcommands.oath.synthesize.SynthCommand());
		addCommand(new gear.subcommands.oath.oathbus.OATHBusCommand());
		addCommand(new gear.subcommands.oath.pick.OathPickCommand());

		addCommand(new gear.subcommands.reml.REMLCommand());
		addCommand(new gear.subcommands.hereg.HERegCommand());
		addCommand(new gear.subcommands.hefam.HEFamCommand());
		addCommand(new gear.subcommands.heregb.HERegBCommand());

		addCommand(new gear.subcommands.fastpca.FastPCACommand());
		
		addCommand(new gear.subcommands.simulationdipop.SimulationDiPopCommand());
		addCommand(new gear.subcommands.simuqtreal.SimulationQTRealCommand());
		addCommand(new gear.subcommands.simulationmpheno.SimulationMPCommand());
		addCommand(new gear.subcommands.simulationfammpheno.SimuFamilyMPCommand());
		addCommand(new gear.subcommands.simulationfam.SimuFamilyCommand());
		addCommand(new gear.subcommands.simulationqt.SimulationQTCommand());
		addCommand(new gear.subcommands.simulationcc.SimulationCCCommand());
		addCommand(new gear.subcommands.simulationinbred.SimulationInbredCommand());

		addCommand(new gear.subcommands.ibd.IBDCommand());		
		addCommand(new gear.subcommands.impute.ImputeCommand());
	}

	private void addCommand(Command cmd) {
		cmdMap.put(cmd.getName(), cmd);
		Iterator<String> aliasIter = cmd.getAliases().iterator();
		while (aliasIter.hasNext()) {
			cmdMap.put(aliasIter.next(), cmd);
		}
	}

	public Command getCommand(String sNameOrAlias) {
		return cmdMap.get(sNameOrAlias);
	}

	public SortedSet<Command> getCommandSet() {
		return new TreeSet<Command>(cmdMap.values());
	}

	private TreeMap<String, Command> cmdMap = new TreeMap<String, Command>();

	public static void main(String[] args) {
		if (args.length == 0) {
			System.out.println("Type 'gear help' or 'java -jar gear.jar help' for usage.");
			System.exit(1);
		}

		subcmdName = args[0];
		String[] subcmdArgs = new String[args.length - 1];
		System.arraycopy(args, 1, subcmdArgs, 0, subcmdArgs.length);

		Command subcmd = INSTANCE.getCommand(subcmdName);

		if (subcmd != null) {
			subcmd.execute(subcmdArgs, subcmdName);
		}
//		} else {
//
//			CmdArgs.INSTANCE.parse(args);
//
//			Logger.setLogFiles(CmdArgs.INSTANCE.out);
//			Logger.hasUserLogTag(false);
//			Logger.printUserLog(AboutInfo.WELCOME_MESSAGE);
//			Logger.hasUserLogTag(true);
//			Logger.printUserLog("Analysis started: " + Calendar.getInstance().getTime() + "\n");
//
//			MonitorThread monitor = new MonitorThread();
//			monitor.start();
//
//			CmdArgs.INSTANCE.printOptionsInEffect();
//
//			if (CmdArgs.INSTANCE.calOption) {
//				H2Transformer H2 = new H2Transformer();
//				H2.H2();
//
//			} else if (CmdArgs.INSTANCE.nontransFlag) {
//				NonTransmitted nt = new NonTransmitted();
//				nt.GenerateNonTransmitted();
//
//			} else if (CmdArgs.INSTANCE.naiveImputFlag) {
//				NaiveImputation ni = new NaiveImputation();
//				ni.Imputation();
//			} else if (CmdArgs.INSTANCE.imputeFlag) {
//				ImputeProbabilityBestGuess bestGuess = new ImputeProbabilityBestGuess();
//				bestGuess.convert();
//			}
//
//			monitor.stopMonitoring();
//
//			Logger.printUserLog("");
//			Logger.printUserLog("Analysis finished: " + Calendar.getInstance().getTime());
//			Logger.printUserLog("Peak memory consumption: " + monitor.getPeakMemoryFormatString());
	}

}
