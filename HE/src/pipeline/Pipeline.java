package pipeline;
import merge.MergeTwoFile;
import write.WriteBedSNPMajor;
import he.H2Transformer;
import he.HERegression;
import realcheck.RealCheck;
import realcheck.RealCheckOne;
import simulation.RealDataSimulation;
import simulation.SimuFamily;
import simulation.SimuPolyCC;
import strand.MakePredictor;
import strand.MakePredictor2;
import strand.Strand;
import sumstat.FrequencyCalculator;
import parameter.Parameter;
import profile.MaCHDosageProfile;
import pscontrol.NonTransmitted;

public class Pipeline {

	public static void main(String[] args) {
		Parameter p = new Parameter();
		p.commandListenor(args);

		if (Parameter.scoreFlag) {
			MaCHDosageProfile mach = new MaCHDosageProfile(p);
			mach.makeProfile();
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
			polyCC.writeFile();
			polyCC.writeLog();

		} else if (Parameter.simupolyQTFlag) {
			
		} else if (Parameter.sumStatFlag) {
			if (Parameter.sumStatOption[Parameter.freq]) {
				FrequencyCalculator fc = new FrequencyCalculator(p);
				fc.CalculateAlleleFrequency();
				System.out.println(fc);

			} else if (Parameter.sumStatOption[Parameter.geno_freq]) {
				FrequencyCalculator fc = new FrequencyCalculator(p);
				fc.CalculateAlleleFrequency();
				System.out.println(fc);				

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

		} else if (Parameter.heFlag){
			HERegression HER = null;
			HER = new HERegression(p);
			HER.Regression();

		}
	}

}
