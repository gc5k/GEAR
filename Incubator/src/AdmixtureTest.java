import java.io.IOException;
import java.util.ArrayList;
import admixture.AdmixtureConstant;
import admixture.ChromosomeGenerator;
import admixture.DNAStirrer;
import admixture.GenerateColony;
import admixture.HotSpot;
import admixture.phenotype.PhenotypeGenerator;
import admixture.phenotype.QualityControl;

/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */

public class AdmixtureTest {

	public static void main(String[] args) {
		// //////// general components
		long seed = 2011;
		int disease_chr = 1;

		// logistic regression

		String[] f = { "0000", "1010", "1111" };
		double[] g_e = { 1, 0.5, 1 };
		int[] chr = { 1, 1 };
		int[] loci = { 0, 1 };
		double mu = 0;
		double dev = Math.sqrt(10);
		PhenotypeGenerator pg = new PhenotypeGenerator(f, g_e, chr, loci, mu, dev);
		HotSpot hs = new HotSpot();

		int N_chr = 2;
		double[] w = {0.8, 0.2};
		ArrayList<DNAStirrer> DNAPool = new ArrayList<DNAStirrer>();
		ArrayList<ChromosomeGenerator> CG = new ArrayList<ChromosomeGenerator>();
		for (int i = 0; i < N_chr; i++) {

			DNAStirrer ds = new DNAStirrer(i, 1, 10000, 5 * (i + 1), AdmixtureConstant.Without_Genetic_Drift, w);
			ds.DNAStir(1);
			DNAPool.add(ds);
			ChromosomeGenerator cg = new ChromosomeGenerator(ds.AdmixedSNPFrequencyPanel());
			cg.setSeed(seed + i);
			CG.add(cg);
		}

		double[] disease_rate = { 0.2, 0.5 };
		int N_phe = 3;
		GenerateColony GC = new GenerateColony(N_phe, seed, disease_chr, disease_rate, hs, DNAPool, CG, pg, AdmixtureConstant.free_recombination);

		// specific components
		// family
		int N_Fam = 100;
		int N_Kid = 2;
		int N_aff_Kid = 1;
		QualityControl qc = new QualityControl(N_aff_Kid, AdmixtureConstant.No_selection);
		GC.GenerateNewFamHab(N_Fam, N_Kid, qc);

		int N_case = 50;
		int N_control = 50;
		QualityControl qc_c = new QualityControl(N_case, N_control, AdmixtureConstant.CaseControl);
		GC.GenerateCCHab(N_case + N_control, 1, qc_c);

		try {
			GC.printGenotype2file("ped.txt", "phe.txt", !AdmixtureConstant.printAllele);
		} catch (IOException e) {
			e.printStackTrace(System.err);
		}
	}
}
