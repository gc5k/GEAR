package admixture.population.scheme;

import java.io.IOException;
import java.util.ArrayList;

import admixture.AdmixtureConstant;
import admixture.population.AlleleFrequencyReader;
import admixture.population.GenerateColony;
import admixture.population.genome.DNAStirrer;
import admixture.population.genome.HotSpot;
import admixture.population.genome.chromosome.ChromosomeGenerator;
import admixture.population.phenotype.PhenotypeGenerator;
import admixture.population.phenotype.QualityControl;


/**
*
* @author Guo-Bo Chen, chenguobo@gmail.com
*/
public class GeneFlowModel {

	/**
	 * @param args
	 */
	public static void main(String[] args) {

		long seed = 2012;
		int control_chr = 0;
		boolean isNullHypothesis = false;
		
		// specific components
		// family
		int N_Fam = 300;
		int N_Kid = 2;
		int N_aff_Kid = 1;
		
		int N_case = 0;
		int N_control = 0;

		// logistic regression
		String[] f = { "0000", "1010", "1111" };
		double[] g_e = { 1, 0.5, 1 };
		int[] chr = { 1, 1 };
		int[] loci = { 0, 1 };
		double mu = 0;
		double dev = Math.sqrt(10);

		double[] disease_rate = { 0.1, 0.3 };
		int N_phe = 3;

		String[] chr_file = {"allele_freq_chr1.txt", "allele_freq_chr2.txt"};
		double[] w = {0.8, 0.2};

		PhenotypeGenerator pg = new PhenotypeGenerator(f, g_e, chr, loci, mu, dev);
		HotSpot hs = new HotSpot();


		ArrayList<DNAStirrer> DNAPool = new ArrayList<DNAStirrer>();
		ArrayList<ChromosomeGenerator> CG = new ArrayList<ChromosomeGenerator>();
		for (int i = 0; i < chr_file.length; i++) {
			AlleleFrequencyReader afr = new AlleleFrequencyReader(chr_file[i]);
			DNAStirrer ds = new DNAStirrer(afr, 1, 10000, AdmixtureConstant.Without_Genetic_Drift, w);
			ds.DNAStir(1);
			DNAPool.add(ds);
			ChromosomeGenerator cg = new ChromosomeGenerator(ds.AncestrySNPFreqencyPanel(), ds.CurrSNPOrigine());
			cg.setSeed(seed + i);
			CG.add(cg);
		}

		GenerateColony GC = new GenerateColony.Builder(disease_rate, hs, DNAPool, CG, pg).numPhenotype(N_phe)
								.seed(seed).diseaseChr(control_chr).recombinationFree(AdmixtureConstant.free_recombination)
								.isNullHypothesis(isNullHypothesis).build();

		QualityControl qc = new QualityControl(N_aff_Kid, AdmixtureConstant.FamilyExactAffected);
		GC.GenerateFamHab(N_Fam, N_Kid, qc);

		QualityControl qc_c = new QualityControl(N_case, N_control, AdmixtureConstant.CaseControl);
		GC.GenerateCCHab(N_case + N_control, 1, qc_c);

		try {
			GC.printGenotype2file("ped.txt", "phe.txt", !AdmixtureConstant.printAllele, 
					!AdmixtureConstant.printLinked);
			GC.printGenotype2file("L_ped.txt", "L_phe.txt", !AdmixtureConstant.printAllele, 
					AdmixtureConstant.printLinked);

			GC.printFamilyGenotype2file("Fam_ped.txt", "Fam_phe.txt", !AdmixtureConstant.printAllele,
					!AdmixtureConstant.printLinked);
			GC.printFamilyGenotype2file("L_Fam_ped.txt", "L_Fam_phe.txt", !AdmixtureConstant.printAllele,
					AdmixtureConstant.printLinked);

			GC.printCCGenotype2file("CC_ped.txt", "CC_phe.txt", !AdmixtureConstant.printAllele,
					!AdmixtureConstant.printLinked);
			GC.printCCGenotype2file("L_CC_ped.txt", "L_CC_phe.txt", !AdmixtureConstant.printAllele,
					AdmixtureConstant.printLinked);
//			GC.printUnrelatedIndividual("PCA_ped.txt", "PCA_phe.txt", !AdmixtureConstant.printAllele, true);
//			GC.printFounder("F_ped.txt", "F_phe.txt", !AdmixtureConstant.printAllele, true);
//			GC.printOffspringCC("KCC_ped.txt", "KCC_phe.txt", !AdmixtureConstant.printAllele);
		} catch (IOException e) {
			e.printStackTrace(System.err);
		}
	}
}
