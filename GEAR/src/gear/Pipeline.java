package gear;

import java.util.Calendar;

import epem.GRMStat;

import merge.MergeTwoFile;
import write.WriteBedSNPMajor;
import gear.profile.MachDosageProfiler;
import gear.profile.ProfilerBase;
import gear.profile.RiskScoreProfiler;
import gear.util.Logger;
import grm.MakeGRM;
import he.HECalculate;
import he.HEPermutation;
import he.HERead;
import he.h2trans.H2Transformer;
import realcheck.RealCheck;
import realcheck.RealCheckOne;
import simulation.RealDataSimulation;
import simulation.SimuFamily;
import simulation.SimuPolyCC;
import simulation.SimuPolyQT;
import strand.MakePredictor;
import strand.MakePredictor2;
import strand.Strand;
import sumstat.FrequencyCalculator;
import sumstat.Inbreeding;
import pscontrol.NonTransmitted;

public class Pipeline
{
	public static void main(String[] args)
	{
		CmdArgs.INSTANCE.commandListener(args);

		Logger.setLogFiles(CmdArgs.INSTANCE.out);
		Logger.printUserLog(AboutInfo.WELCOME_MESSAGE);
		Logger.printUserLog("Analysis started: " + Calendar.getInstance().getTime() + "\n");
		
		// Print the options set.
		Logger.printUserLog("Options received: ");
		StringBuffer sb = null;
		for (int i = 0; i < args.length; i++)
		{
			if (args[i].startsWith("--"))
			{
				if (i > 1)
				{
					Logger.printUserLog(sb.toString());
				}
				sb = new StringBuffer("  ");
			}
			sb.append(args[i] + " ");
		}
		Logger.printUserLog(sb.toString());
		Logger.printUserLog("");

		if (CmdArgs.INSTANCE.getHpcArgs().isSet())
		{
			HPC.genScript(args);
		}
		else if (CmdArgs.INSTANCE.getProfileArgs().isSet())
		{
			ProfilerBase profiler;
			if (CmdArgs.INSTANCE.getFileArgs().isSet() || CmdArgs.INSTANCE.getBFileArgs(0).isSet())
			{
				Logger.printUserLog("Generating risk profiles for genotypes.");
				profiler = new RiskScoreProfiler();
			}
			else
			{
				Logger.printUserLog("Generating risk profiles for MaCH dosage.");
				profiler = new MachDosageProfiler();
			}
			profiler.makeProfile();
		}
		else if (CmdArgs.INSTANCE.hasStrandOption())
		{
			Strand strand = new Strand();
			strand.Merge();
		}
		else if (CmdArgs.INSTANCE.hasMergeOption())
		{
			MergeTwoFile mtf = new MergeTwoFile();
			mtf.Merge();
		}
		else if (CmdArgs.INSTANCE.hasMakePredictorOption())
		{
			MakePredictor mp = new MakePredictor();
			mp.BuildPredictor();
		}
		else if (CmdArgs.INSTANCE.hasMakePredictor2Option())
		{
			MakePredictor2 mp2 = new MakePredictor2();
			mp2.BuildPredictor();
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
			SimuFamily simuFam = new SimuFamily(
					CmdArgs.INSTANCE.simu_fam_size,
					CmdArgs.INSTANCE.simu_fam_marker, CmdArgs.INSTANCE.seed);
			simuFam.generateSample();

		}
		else if (CmdArgs.INSTANCE.simuRealData)
		{
			RealDataSimulation rdSimu = new RealDataSimulation();
			rdSimu.GenerateSample();

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
				Logger.printUserLog(fc.toString());

			} else if (CmdArgs.INSTANCE.genoFreqFlag)
			{
				FrequencyCalculator fc = new FrequencyCalculator();
				fc.CalculateAlleleFrequency();
				Logger.printUserLog(fc.toString());

			} else if (CmdArgs.INSTANCE.fstFlag)
			{
				Inbreeding inb = new Inbreeding();
				inb.CalculateFst();
			}
		}
		else if (CmdArgs.INSTANCE.makebedFlag)
		{
			WriteBedSNPMajor bedWriter = new WriteBedSNPMajor();
			bedWriter.WriteFile();

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
			HERead hr = new HERead();
			HECalculate HC = new HECalculate(hr);
			HC.Regression();

			if (CmdArgs.INSTANCE.permFlag)
			{
				HEPermutation hp = new HEPermutation(hr);
				hp.Permutation();
			}
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
		Logger.printUserLog("");
		Logger.printUserLog("Analysis finished: "
				+ Calendar.getInstance().getTime());
	}

}
