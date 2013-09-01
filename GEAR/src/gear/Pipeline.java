package gear;

import java.util.Calendar;


import gear.ancestry.GenPCA;
import gear.ancestry.HapMapPCA;
import gear.epem.GRMStat;
import gear.grm.MakeGRM;
import gear.he.HEMCalculate;
import gear.he.HEMRead;
import gear.he.h2trans.H2Transformer;
import gear.imputation.NaiveImputation;
import gear.profile.Profiler;
import gear.pscontrol.NonTransmitted;
import gear.realcheck.RealCheck;
import gear.realcheck.RealCheckOne;
import gear.simulation.RealDataSimulation;
import gear.simulation.RealDataSimulationQT;
import gear.simulation.SimuFamily;
import gear.simulation.SimuPolyCC;
import gear.simulation.SimuPolyQT;
import gear.sumstat.FrequencyCalculator;
import gear.sumstat.Inbreeding;
import gear.util.Logger;
import gear.util.MonitorThread;
import gear.write.WriteBedSNPMajor;

public class Pipeline
{
	public static void main(String[] args)
	{
		if (args.length == 0)
		{
			System.out.println("Type 'gear help' or 'java -jar gear.jar help' for usage.");
			System.exit(1);
		}
		
		CmdArgs.INSTANCE.parse(args);
		
		Logger.setLogFiles(CmdArgs.INSTANCE.out);
		Logger.printUserLog(AboutInfo.WELCOME_MESSAGE);
		Logger.printUserLog("Analysis started: " + Calendar.getInstance().getTime() + "\n");
		
		MonitorThread monitor = new MonitorThread();
		monitor.start();

		CmdArgs.INSTANCE.printOptionsInEffect();

		if (CmdArgs.INSTANCE.getHpcArgs().isSet())
		{
			HPC.genScript(args);
		}
		else if (CmdArgs.INSTANCE.getProfileArgs().isSet())
		{
			Profiler.makeProfile();
		}
		else if (CmdArgs.INSTANCE.hasRealCheckOption())
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
		else if (CmdArgs.INSTANCE.simufamFlag)
		{
			SimuFamily simuFam = new SimuFamily();
			simuFam.generateSample();

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
		else if (CmdArgs.INSTANCE.pcaFlag)
		{
			GenPCA gpca = new GenPCA();
		}
		else if (CmdArgs.INSTANCE.hppcaFlag)
		{
			HapMapPCA hppca = new HapMapPCA();
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
		
		monitor.stopMonitoring();
		
		Logger.printUserLog("");
		Logger.printUserLog("Analysis finished: " + Calendar.getInstance().getTime());
		Logger.printUserLog("Peak memory consumption: " + monitor.getPeakMemoryFormatString());
	}
}
