package pipeline;
import java.util.Calendar;

import merge.MergeTwoFile;
import write.WriteBedSNPMajor;
import he.HECalculate;
import he.HEPermutation;
import he.HERead;
import he.HERegression;
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
import parameter.Parameter;
import profile.MaCHDosageProfile;
import profile.RiskScore;
import pscontrol.NonTransmitted;

public class Pipeline {

	public static void main(String[] args) {
		Parameter.INSTANCE.commandListener(args);

		System.out.print(Parameter.version);
		Calendar calendar = Calendar.getInstance();
		System.out.println("\nThe analysis started at: " + calendar.getTime() + "\n");

		if (Parameter.INSTANCE.scoreFlag) {
			if (Parameter.INSTANCE.hasBFileOption()) {
				RiskScore rs = new RiskScore();
				rs.makeProfile();
			} else {
				MaCHDosageProfile mach = new MaCHDosageProfile();
				mach.makeProfile();
			}
		} else if (Parameter.INSTANCE.hasStrandOption()) {
			Strand strand = new Strand();
			strand.Merge();
		} else if (Parameter.INSTANCE.hasMergeOption()) {
			MergeTwoFile mtf = new MergeTwoFile();
			mtf.Merge();
		} else if (Parameter.INSTANCE.hasMakePredictorOption()) {
			MakePredictor mp = new MakePredictor();
			mp.BuildPredictor();
		} else if (Parameter.INSTANCE.hasMakePredictor2Option()) {
			MakePredictor2 mp2 = new MakePredictor2();
			mp2.BuildPredictor();
		} else if (Parameter.INSTANCE.hasRealCheckOption()) {
			if (Parameter.INSTANCE.hasBFile2Option()) {
				RealCheck realcheck = new RealCheck();
				realcheck.Check();
			} else {
				RealCheckOne realcheckone = new RealCheckOne();
				realcheckone.Check();
			}
		} else if (Parameter.simufamFlag) {
			SimuFamily simuFam = new SimuFamily (Parameter.INSTANCE.simu_fam_size,
												 Parameter.INSTANCE.simu_fam_marker,
												 Parameter.INSTANCE.seed);
			simuFam.generateSample();

		} else if (Parameter.simuRealData) {
			RealDataSimulation rdSimu = new RealDataSimulation();
			rdSimu.GenerateSample();
			
		}  else if (Parameter.simupolyCCFlag) {
			SimuPolyCC polyCC = new SimuPolyCC();
			polyCC.GenerateSample();

		} else if (Parameter.simupolyQTFlag) {
			SimuPolyQT polyQT = new SimuPolyQT();
			polyQT.generateSample();
			
		} else if (Parameter.sumStatFlag) {
			if (Parameter.freqFlag) {
				FrequencyCalculator fc = new FrequencyCalculator();
				fc.CalculateAlleleFrequency();
				System.out.println(fc);

			} else if (Parameter.genoFreqFlag) {
				FrequencyCalculator fc = new FrequencyCalculator();
				fc.CalculateAlleleFrequency();
				System.out.println(fc);				

			} else if (Parameter.fstFlag ) {
				Inbreeding inb = new Inbreeding();
				inb.CalculateFst();
			}
		} else if (Parameter.makebedFlag) {
			WriteBedSNPMajor bedWriter = new WriteBedSNPMajor();
			bedWriter.WriteFile();

		} else if (Parameter.calOption) {
			H2Transformer H2 = null;
			H2 = new H2Transformer();
			H2.H2();

		} else if (Parameter.nontransFlag) {
			NonTransmitted nt = new NonTransmitted();
			nt.GenerateNonTransmitted();

		} else if (Parameter.heFlag) {
			HERead hr = new HERead();
			HECalculate HC = new HECalculate(hr);
			HC.Regression();

			if(Parameter.INSTANCE.permFlag) {
				HEPermutation hp = new HEPermutation(hr);
				hp.Permutation();
			}
//			HERegression HER = null;
//			HER = new HERegression(p);
//			HER.Regression();

		}
		
		System.out.println("\nThe analysis ended at: " + calendar.getTime() + "\n");

	}

}
