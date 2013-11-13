package gear;

import java.util.Calendar;
import java.util.Iterator;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import gear.epem.GRMStat;
import gear.grm.MakeGRM;
import gear.he.HEMCalculate;
import gear.he.HEMRead;
import gear.he.h2trans.H2Transformer;
import gear.imputation.NaiveImputation;
import gear.metawatchdog.MetaWatchdog;
import gear.pscontrol.NonTransmitted;
import gear.realcheck.RealCheck;
import gear.realcheck.RealCheckOne;
import gear.simulation.RealDataSimulation;
import gear.simulation.RealDataSimulationQT;
import gear.simulation.SimuPolyCC;
import gear.simulation.SimuPolyQT;
import gear.sumstat.FrequencyCalculator;
import gear.sumstat.Inbreeding;
import gear.util.Logger;
import gear.util.MonitorThread;
import gear.write.WriteBedSNPMajor;

public enum Gear
{
	INSTANCE;
	
	private Gear()
	{
		addCommand(new gear.bluppca.BlupPcaCommand());
		addCommand(new gear.encrypt.EnigmaCommand());
		addCommand(new gear.help.HelpCommand());
		addCommand(new gear.hpc.HpcCommand());
		addCommand(new gear.profile.ProfileCommand());
		addCommand(new gear.simulation.SimuFamilyCommand());
	}
	
	private void addCommand(Command cmd)
	{
		cmdMap.put(cmd.getName(), cmd);
		Iterator<String> aliasIter = cmd.getAliases().iterator();
		while (aliasIter.hasNext())
		{
			cmdMap.put(aliasIter.next(), cmd);
		}
	}
	
	public Command getCommand(String sNameOrAlias)
	{
		return cmdMap.get(sNameOrAlias);
	}
	
	public SortedSet<Command> getCommandSet()
	{
		 return new TreeSet<Command>(cmdMap.values());
	}
	
	private TreeMap<String, Command> cmdMap = new TreeMap<String, Command>();
	
	public static void main(String[] args)
	{
		if (args.length == 0)
		{
			System.out.println("Type 'gear help' or 'java -jar gear.jar help' for usage.");
			System.exit(1);
		}
		
		String subcmdName = args[0];
		String[] subcmdArgs = new String[args.length - 1];
		System.arraycopy(args, 1, subcmdArgs, 0, subcmdArgs.length);
		
		Command subcmd = INSTANCE.getCommand(subcmdName);
		
		if (subcmd != null)
		{
			subcmd.execute(subcmdArgs);
		}
		else
		{
			// TODO: Move the following routines to the subcommand system and print 'Unknown subcommand'.
			
			CmdArgs.INSTANCE.parse(args);
			
			Logger.setLogFiles(CmdArgs.INSTANCE.out);
			Logger.hasUserLogTag(false);
			Logger.printUserLog(AboutInfo.WELCOME_MESSAGE);
			Logger.hasUserLogTag(true);
			Logger.printUserLog("Analysis started: " + Calendar.getInstance().getTime() + "\n");
			
			MonitorThread monitor = new MonitorThread();
			monitor.start();
	
			CmdArgs.INSTANCE.printOptionsInEffect();
	
			if (CmdArgs.INSTANCE.hasRealCheckOption())
			{
				if (CmdArgs.INSTANCE.getBFileArgs(1).isSet())
				{
					RealCheck realcheck = new RealCheck();
					realcheck.Check();
				}
				else
				{
					RealCheckOne realcheckone = new RealCheckOne();
					realcheckone.Check();
				}
			}
			else if (CmdArgs.INSTANCE.bsimuFlag)
			{
				if (CmdArgs.INSTANCE.simupolyCCFlag) 
				{
					RealDataSimulation rdSimu = new RealDataSimulation();
					rdSimu.GenerateSample();				
				}
				else 
				{
					RealDataSimulationQT rdSimuQT = new RealDataSimulationQT();
					rdSimuQT.GenerateSample();
				}
			}
			else if (CmdArgs.INSTANCE.simupolyCCFlag)
			{
				SimuPolyCC polyCC = new SimuPolyCC();
				polyCC.GenerateSample();
	
			}
			else if (CmdArgs.INSTANCE.simupolyQTFlag)
			{
				SimuPolyQT polyQT = new SimuPolyQT();
				polyQT.generateSample();
	
			}
			else if (CmdArgs.INSTANCE.sumStatFlag)
			{
				if (CmdArgs.INSTANCE.freqFlag)
				{
					FrequencyCalculator fc = new FrequencyCalculator();
					fc.CalculateAlleleFrequency();
					fc.PrintOut();
	
				} else if (CmdArgs.INSTANCE.genoFreqFlag)
				{
					FrequencyCalculator fc = new FrequencyCalculator();
					fc.CalculateAlleleFrequency();
					fc.PrintOut();
	
				} else if (CmdArgs.INSTANCE.fstFlag)
				{
					Inbreeding inb = new Inbreeding();
					inb.CalculateFst();
				}
			}
			else if (CmdArgs.INSTANCE.calOption)
			{
				H2Transformer H2 = new H2Transformer();
				H2.H2();
	
			}
			else if (CmdArgs.INSTANCE.nontransFlag)
			{
				NonTransmitted nt = new NonTransmitted();
				nt.GenerateNonTransmitted();
	
			}
			else if (CmdArgs.INSTANCE.hasHEOption())
			{
				HEMRead mhr = new HEMRead();
				HEMCalculate mhc = new HEMCalculate(mhr);
				mhc.Regression();
	//			if (CmdArgs.INSTANCE.getHEArgs().isSingleGrm())
	//			{
	//				HERead hr = new HERead();
	//				HECalculate hc = new HECalculate(hr);
	//				hc.Regression();
	//
	//				if (CmdArgs.INSTANCE.permFlag)
	//				{
	//					HEPermutation hp = new HEPermutation(hr);
	//					hp.Permutation();
	//				}
	//			}
			}
			else if (CmdArgs.INSTANCE.grmstatFlag)
			{
				GRMStat gs = new GRMStat();
				gs.GetGRMStats();
	
			}
			else if (CmdArgs.INSTANCE.GRMFlag)
			{
				MakeGRM mg = new MakeGRM();
				if (CmdArgs.INSTANCE.grmRangeFlag)
				{
					mg.makeGeneticRelationshipScore(
							CmdArgs.INSTANCE.grm_range[0],
							CmdArgs.INSTANCE.grm_range[1]);
				}
				else if (CmdArgs.INSTANCE.grmPartitionFlag)
				{
					mg.GRMPartitioning(args);
				}
				else
				{
					mg.makeGeneticRelationshipScore();
				}
			}
			else if (CmdArgs.INSTANCE.naiveImputFlag)
			{
				NaiveImputation ni = new NaiveImputation();
				ni.Imputation();
			}
			else if (CmdArgs.INSTANCE.makebedFlag)
			{
				WriteBedSNPMajor bedWriter = new WriteBedSNPMajor();
				bedWriter.WriteFile();
	
			}
			else if (CmdArgs.INSTANCE.watchdogFlag)
			{
				MetaWatchdog wd = new MetaWatchdog();
				wd.Bark();
			}
			
			monitor.stopMonitoring();
			
			Logger.printUserLog("");
			Logger.printUserLog("Analysis finished: " + Calendar.getInstance().getTime());
			Logger.printUserLog("Peak memory consumption: " + monitor.getPeakMemoryFormatString());
		}
	}
}
