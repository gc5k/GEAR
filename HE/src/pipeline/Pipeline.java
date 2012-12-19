package pipeline;
import java.util.Calendar;

import merge.MergeTwoFile;
import write.WriteBedSNPMajor;
import he.H2Transformer;
import he.HECalculate;
import he.HEPermutation;
import he.HERead;
import he.HERegression;
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
		Parameter p = new Parameter();
		p.commandListenor(args);

		System.out.print(Parameter.version);
		Calendar calendar = Calendar.getInstance();
		System.out.println("\nThe analysis started at: " + calendar.getTime() + "\n");

		if (Parameter.scoreFlag) {
			if (Parameter.bfileOption) {
				RiskScore rs = new RiskScore(p);
				rs.makeProfile();
			} else {
				MaCHDosageProfile mach = new MaCHDosageProfile(p);
				mach.makeProfile();
			}
		} else if (Parameter.strandFlag) {
			Strand strand = new Strand(p);
			strand.Merge();
		} else if (Parameter.mergeFlag) {
			MergeTwoFile mtf = new MergeTwoFile(p);
			mtf.Merge();
		} else if (Parameter.makePredictorFlag) {
			MakePredictor mp = new MakePredictor(p);
			mp.BuildPredictor();
		} else if (Parameter.makePredictor2Flag) {
			MakePredictor2 mp2 = new MakePredictor2(p);
			mp2.BuildPredictor();
		} else if (Parameter.realcheckFlag) {
			if(Parameter.bfile2 != null) {
				RealCheck realcheck = new RealCheck(p);
				realcheck.Check();
			} else {
				RealCheckOne realcheckone = new RealCheckOne(p);
				realcheckone.Check();
			}
		} else if (Parameter.simufamFlag) {
			SimuFamily simuFam = new SimuFamily(p);
			simuFam.generateSample();

		} else if (Parameter.simuRealData) {
			RealDataSimulation rdSimu = new RealDataSimulation(p);
			rdSimu.GenerateSample();
			
		}  else if (Parameter.simupolyCCFlag) {
			SimuPolyCC polyCC = new SimuPolyCC(p);
			polyCC.GenerateSample();

		} else if (Parameter.simupolyQTFlag) {
			SimuPolyQT polyQT = new SimuPolyQT(p);
			polyQT.generateSample();
			
		} else if (Parameter.sumStatFlag) {
			if (Parameter.freqFlag) {
				FrequencyCalculator fc = new FrequencyCalculator(p);
				fc.CalculateAlleleFrequency();
				System.out.println(fc);

			} else if (Parameter.genoFreqFlag) {
				FrequencyCalculator fc = new FrequencyCalculator(p);
				fc.CalculateAlleleFrequency();
				System.out.println(fc);				

			} else if (Parameter.fstFlag ) {
				Inbreeding inb = new Inbreeding(p);
				inb.CalculateFst();
			}
		} else if (Parameter.makebedFlag) {
			WriteBedSNPMajor bedWriter = new WriteBedSNPMajor(p);
			bedWriter.WriteFile();

		} else if (Parameter.calOption) {
			H2Transformer H2 = null;
			H2 = new H2Transformer(p);
			H2.H2();

		} else if (Parameter.nontransFlag) {
			NonTransmitted nt = new NonTransmitted(p);
			nt.GenerateNonTransmitted();

		} else if (Parameter.heFlag) {
			HERead hr = new HERead(p);
			HECalculate HC = new HECalculate(hr);
			HC.Regression();

			if(p.permFlag) {
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
