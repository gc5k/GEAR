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
		int N_Fam = 50;
		int N_Kid = 2;
		int N_aff_Kid = 1;

		int N_case = 50;
		int N_control = 50;

		// logistic regression
		String[] f = { "0000", "1010", "1111" };
		double[] g_e = { 1, 0.5, 1 };
		int[] chr = { 1, 1 };
		int[] loci = { 0, 1 };
		double mu = 0;
		double dev = Math.sqrt(10);

		double[] disease_rate = { 0, 0 };
		int N_phe = 3;

		String[] chr_file = { "allele_freq_chr1.txt", "allele_freq_chr2.txt" };
		int[] AIM_number = {0, 0};
		double[] w = { 0.8, 0.2 };

		PhenotypeGenerator pg = new PhenotypeGenerator(f, g_e, chr, loci, mu, dev);
		HotSpot hs = new HotSpot();

		ArrayList<DNAStirrer> DNAPool = new ArrayList<DNAStirrer>();
		ArrayList<ChromosomeGenerator> CG = new ArrayList<ChromosomeGenerator>();
		for (int i = 0; i < chr_file.length; i++) {
			AlleleFrequencyReader afr = new AlleleFrequencyReader(chr_file[i], AIM_number[i]);
			DNAStirrer ds = new DNAStirrer(afr, 1, 10000, AdmixtureConstant.Without_Genetic_Drift, w);
			ds.DNAStir(1);
			DNAPool.add(ds);
			ChromosomeGenerator cg = new ChromosomeGenerator(ds.AncestrySNPFreqencyPanel(), ds.CurrSNPOrigine());
			cg.setSeed(seed + i);
			CG.add(cg);
		}

		GenerateColony GC = new GenerateColony.Builder().diesaseRate(disease_rate).hotSpot(hs).numPhenotype(N_phe)
				.ChrGenerator(CG).DNAPool(DNAPool).PheGenerator(pg).seed(seed).diseaseChr(control_chr)
				.isNullHypothesis(isNullHypothesis).build();

		for (int rep = 0; rep < 2; rep++) {
			QualityControl qc = new QualityControl(N_aff_Kid, AdmixtureConstant.FamilyExactAffected);
			QualityControl qc_c = new QualityControl(N_case, N_control, AdmixtureConstant.CaseControl);
			
			GenerateColony.setCurrFamilyID(0);
			GC.GenerateFamHab(N_Fam, N_Kid, qc);
			GC.GenerateCCHab(N_case + N_control, 1, qc_c);

			try {
				StringBuilder Gsb = new StringBuilder("ped.txt");
				Gsb.insert(0, rep);
				StringBuilder Psb = new StringBuilder("phe.txt");
				Psb.insert(0, rep);
				StringBuilder L_Gsb = new StringBuilder("L_ped.txt");
				L_Gsb.insert(0, rep);
				StringBuilder L_Psb = new StringBuilder("L_phe.txt");
				L_Psb.insert(0, rep);
				GC.printGenotype2file(Gsb.toString(), Psb.toString(), !AdmixtureConstant.printAllele,
						!AdmixtureConstant.printLinked);
				GC.printGenotype2file(L_Gsb.toString(), L_Psb.toString(), AdmixtureConstant.printAllele,
						AdmixtureConstant.printLinked);

				StringBuilder F_Gsb = new StringBuilder("Fam_ped.txt");
				F_Gsb.insert(0, rep);
				StringBuilder F_Psb = new StringBuilder("Fam_phe.txt");
				F_Psb.insert(0, rep);
				StringBuilder L_F_Gsb = new StringBuilder("L_Fam_ped.txt");
				L_F_Gsb.insert(0, rep);
				StringBuilder L_F_Psb = new StringBuilder("L_Fam_phe.txt");
				L_F_Psb.insert(0, rep);
				GC.printFamilyGenotype2file(F_Gsb.toString(), F_Psb.toString(), !AdmixtureConstant.printAllele,
						!AdmixtureConstant.printLinked);
				GC.printFamilyGenotype2file(L_F_Gsb.toString(), L_F_Psb.toString(), AdmixtureConstant.printAllele,
						AdmixtureConstant.printLinked);

				StringBuilder C_Gsb = new StringBuilder("CC_ped.txt");
				C_Gsb.insert(0, rep);
				StringBuilder C_Psb = new StringBuilder("CC_phe.txt");
				C_Psb.insert(0, rep);
				StringBuilder L_C_Gsb = new StringBuilder("L_CC_ped.txt");
				L_C_Gsb.insert(0, rep);
				StringBuilder L_C_Psb = new StringBuilder("L_CC_phe.txt");
				L_C_Psb.insert(0, rep);
				GC.printCCGenotype2file(C_Gsb.toString(), C_Psb.toString(), !AdmixtureConstant.printAllele,
						!AdmixtureConstant.printLinked);
				GC.printCCGenotype2file(L_C_Gsb.toString(), L_C_Psb.toString(), AdmixtureConstant.printAllele,
						AdmixtureConstant.printLinked);

				// GC.printUnrelatedIndividual("PCA_ped.txt", "PCA_phe.txt", !AdmixtureConstant.printAllele, true);
				// GC.printFounder("F_ped.txt", "F_phe.txt", !AdmixtureConstant.printAllele, true);
				// GC.printOffspringCC("KCC_ped.txt", "KCC_phe.txt", !AdmixtureConstant.printAllele);
			} catch (IOException e) {
				e.printStackTrace(System.err);
			}
		}
	}
}
