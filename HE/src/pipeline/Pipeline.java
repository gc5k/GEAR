package pipeline;
import write.WriteBedSNPMajor;
import he.H2Transformer;
import he.HERegression;
import realcheck.RealCheck;
import simulation.RealDataSimulation;
import simulation.SimuFamily;
import simulation.SimuPolyCC;
import sumstat.FrequencyCalculator;
import parameter.Parameter;
import pscontrol.NonTransmitted;

public class Pipeline {

	public static void main(String[] args) {
		Parameter p = new Parameter();
		p.commandListenor(args);

		if (Parameter.realcheckFlag) {
			RealCheck realcheck = new RealCheck(p);
			realcheck.Check();
			
		} else if (Parameter.simufamFlag) {
			SimuFamily simuFam = new SimuFamily(p);
			simuFam.generateSample();
		} else if (Parameter.simuFlag) {
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

		} else {
			HERegression HER = null;
			HER = new HERegression(p);
			HER.Regression();

		}
	}

}
