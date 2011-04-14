package admixture.population.scheme;

import java.io.IOException;
import java.util.ArrayList;

import admixture.AdmixtureConstant;
import admixture.population.AlleleFrequencyReader;
import admixture.population.GenerateColony;
import admixture.population.Print2File;
import admixture.population.genome.DNAStirrer;
import admixture.population.genome.HotSpot;
import admixture.population.genome.chromosome.ChromosomeGenerator;
import admixture.population.phenotype.PhenotypeGenerator;
import admixture.population.phenotype.QualityControl;

public class DiscreteAncestryModel {

	public static void main(String[] args) {

		long seed = 2011;
		int control_chr = 0;
		boolean isNullHypothesis = true;
		// logistic regression

		String[] f = { "1111", "2121", "2222" };
		double[] g_e = { 1, 0.5, 1 };
		int[] chr = { 1, 1 };
		int[] loci = { 0, 1 };
		double mu = 0;
		double dev = Math.sqrt(10);
		PhenotypeGenerator pg = new PhenotypeGenerator(f, g_e, chr, loci, mu, dev);
		HotSpot hs = new HotSpot();

		int N_phe = 2;

		// specify AA parameters
		double[] AA_disease_rate = { 0.2 };
		int AA_N_Fam = 100;
		int AA_N_Kid = 2;
		int AA_N_aff_Kid = 1;

		int AA_N_case = 50;
		int AA_N_control = 50;

		String[] AA_chr_file = { "AIM_AA_chr1.txt", "allele_freq_AA_chr2_simu.txt" };
		int[] AIM_number = {0, 0};
		double[] AA_w = { 1 };
		ArrayList<DNAStirrer> AA_DNAPool = new ArrayList<DNAStirrer>();
		ArrayList<ChromosomeGenerator> AA_CG = new ArrayList<ChromosomeGenerator>();
		for (int i = 0; i < AA_chr_file.length; i++) {
			AlleleFrequencyReader afr = new AlleleFrequencyReader(AA_chr_file[i], AIM_number[i]);
			DNAStirrer ds = new DNAStirrer(afr, 1, 10000, AdmixtureConstant.Without_Genetic_Drift, AA_w);
			ds.DNAStir(1);
			AA_DNAPool.add(ds);
			ChromosomeGenerator cg = new ChromosomeGenerator(ds.AncestrySNPFreqencyPanel(), ds.CurrSNPOrigine());
			cg.setSeed(seed + i);
			AA_CG.add(cg);
		}
		GenerateColony AA_GenColony = new GenerateColony(N_phe, seed, control_chr, AA_disease_rate, hs, AA_DNAPool,
				AA_CG, pg, AdmixtureConstant.free_recombination, isNullHypothesis);

		// specify EA parameters
		double[] EA_disease_rate = { 0.5 };
		int EA_N_Fam = 100;
		int EA_N_Kid = 2;
		int EA_N_aff_Kid = 1;

		int EA_N_case = 50;
		int EA_N_control = 50;

		// generate EA
		String[] EA_chr_file = { "AIM_EA_chr1.txt", "allele_freq_EA_chr2_simu.txt" };
		double[] w = { 1 };
		ArrayList<DNAStirrer> EA_DNAPool = new ArrayList<DNAStirrer>();
		ArrayList<ChromosomeGenerator> EA_CG = new ArrayList<ChromosomeGenerator>();
		for (int i = 0; i < EA_chr_file.length; i++) {
			AlleleFrequencyReader afr = new AlleleFrequencyReader(EA_chr_file[i], AIM_number[i]);
			DNAStirrer ds = new DNAStirrer(afr, 1, 10000, AdmixtureConstant.Without_Genetic_Drift, w);
			ds.DNAStir(1);
			EA_DNAPool.add(ds);
			ChromosomeGenerator cg = new ChromosomeGenerator(ds.AncestrySNPFreqencyPanel(), ds.CurrSNPOrigine());
			cg.setSeed(seed + i);
			EA_CG.add(cg);
		}

		boolean printLinked = true;
		for (int rep = 0; rep < 10; rep++) {
			QualityControl AA_qc = new QualityControl(AA_N_aff_Kid, AdmixtureConstant.FamilyExactAffected);
			QualityControl AA_qc_c = new QualityControl(AA_N_case, AA_N_control, AdmixtureConstant.CaseControl);
			
			QualityControl EA_qc = new QualityControl(EA_N_aff_Kid, AdmixtureConstant.FamilyExactAffected);
			QualityControl EA_qc_c = new QualityControl(EA_N_case, EA_N_control, AdmixtureConstant.CaseControl);

			GenerateColony EA_GenColony = new GenerateColony(N_phe, seed, control_chr, EA_disease_rate, hs, EA_DNAPool,
					EA_CG, pg, AdmixtureConstant.free_recombination, isNullHypothesis);

			GenerateColony.setCurrFamilyID(0);
			AA_GenColony.GenerateFamHab(AA_N_Fam, AA_N_Kid, AA_qc);
			AA_GenColony.GenerateCCHab(AA_N_case + AA_N_control, 1, AA_qc_c);

			EA_GenColony.GenerateFamHab(EA_N_Fam, EA_N_Kid, EA_qc);
			EA_GenColony.GenerateCCHab(EA_N_case + EA_N_control, 1, EA_qc_c);

			Print2File p2f = new Print2File();
			p2f.addColony(AA_GenColony);
			p2f.addColony(EA_GenColony);
			try {
				StringBuilder Gsb = new StringBuilder("ped.txt");
				Gsb.insert(0, rep);
				StringBuilder Psb = new StringBuilder("phe.txt");
				Psb.insert(0, rep);
				StringBuilder L_Gsb = new StringBuilder("L_ped.txt");
				L_Gsb.insert(0, rep);
				StringBuilder L_Psb = new StringBuilder("L_phe.txt");
				L_Psb.insert(0, rep);
				p2f.printGenotype2file(Gsb.toString(), Psb.toString(), !AdmixtureConstant.printAllele,
						!printLinked);
				p2f.printGenotype2file(L_Gsb.toString(), L_Psb.toString(), AdmixtureConstant.printAllele,
						printLinked);

				StringBuilder F_Gsb = new StringBuilder("Fam_ped.txt");
				F_Gsb.insert(0, rep);
				StringBuilder F_Psb = new StringBuilder("Fam_phe.txt");
				F_Psb.insert(0, rep);
				StringBuilder L_F_Gsb = new StringBuilder("L_Fam_ped.txt");
				L_F_Gsb.insert(0, rep);
				StringBuilder L_F_Psb = new StringBuilder("L_Fam_phe.txt");
				L_F_Psb.insert(0, rep);
				p2f.printFamilyGenotype2file(F_Gsb.toString(), F_Psb.toString(), !AdmixtureConstant.printAllele,
						!printLinked);
				p2f.printFamilyGenotype2file(L_F_Gsb.toString(), L_F_Psb.toString(), AdmixtureConstant.printAllele,
						printLinked);

				StringBuilder C_Gsb = new StringBuilder("CC_ped.txt");
				C_Gsb.insert(0, rep);
				StringBuilder C_Psb = new StringBuilder("CC_phe.txt");
				C_Psb.insert(0, rep);
				StringBuilder L_C_Gsb = new StringBuilder("L_CC_ped.txt");
				L_C_Gsb.insert(0, rep);
				StringBuilder L_C_Psb = new StringBuilder("L_CC_phe.txt");
				L_C_Psb.insert(0, rep);
				p2f.printCCGenotype2file(C_Gsb.toString(), C_Psb.toString(), !AdmixtureConstant.printAllele,
						!printLinked);
				p2f.printCCGenotype2file(L_C_Gsb.toString(), L_C_Psb.toString(), AdmixtureConstant.printAllele,
						printLinked);

				// p2f.printUnrelatedIndividual("PCA_ped.txt", "PCA_phe.txt", !AdmixtureConstant.printAllele, true);
				// p2f.printFounder("F_ped.txt", "F_phe.txt", !AdmixtureConstant.printAllele, true);
			} catch (IOException e) {
				e.printStackTrace(System.err);
			}
		}
	}
}
