package admixture;

import java.io.IOException;
import java.util.ArrayList;

import admixture.phenotype.PhenotypeGenerator;
import admixture.phenotype.QualityControl;
import admixture.print.Print2File;

public class DiscreteAncestryModel {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
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
		
		String[] AA_chr_file = {"allele_freq_AA_chr1_simu.txt", "allele_freq_AA_chr2_simu.txt"};
		double[] AA_w = {1};
		ArrayList<DNAStirrer> AA_DNAPool = new ArrayList<DNAStirrer>();
		ArrayList<ChromosomeGenerator> AA_CG = new ArrayList<ChromosomeGenerator>();
		for (int i = 0; i < AA_chr_file.length; i++) {
			AlleleFrequencyReader afr = new AlleleFrequencyReader(AA_chr_file[i]);
			DNAStirrer ds = new DNAStirrer(afr, 1, 10000, AdmixtureConstant.Without_Genetic_Drift, AA_w);
			ds.DNAStir(1);
			AA_DNAPool.add(ds);
			ChromosomeGenerator cg = new ChromosomeGenerator(ds.AncestrySNPFreqencyPanel(), ds.CurrSNPOrigine());
			cg.setSeed(seed + i);
			AA_CG.add(cg);
		}

		int N_phe = 3;
		
		double[] AA_disease_rate = { 0.2};
		GenerateColony AA_GC = new GenerateColony(N_phe, seed, disease_chr, AA_disease_rate, 
				hs, AA_DNAPool, AA_CG, pg, AdmixtureConstant.free_recombination);

		// specific components
		// family
		int AA_N_Fam = 100;
		int AA_N_Kid = 2;
		int AA_N_aff_Kid = 1;
		QualityControl AA_qc = new QualityControl(AA_N_aff_Kid, AdmixtureConstant.No_selection);
		AA_GC.GenerateNewFamHab(AA_N_Fam, AA_N_Kid, AA_qc);

		int AA_N_case = 50;
		int AA_N_control = 50;
		QualityControl AA_qc_c = new QualityControl(AA_N_case, AA_N_control, AdmixtureConstant.No_selection);
		AA_GC.GenerateCCHab(AA_N_case + AA_N_control, 1, AA_qc_c);

//generate EA		
		String[] EA_chr_file = {"allele_freq_EA_chr1_simu.txt", "allele_freq_EA_chr2_simu.txt"};
		double[] w = {1};
		ArrayList<DNAStirrer> EA_DNAPool = new ArrayList<DNAStirrer>();
		ArrayList<ChromosomeGenerator> EA_CG = new ArrayList<ChromosomeGenerator>();
		for (int i = 0; i < EA_chr_file.length; i++) {
			AlleleFrequencyReader afr = new AlleleFrequencyReader(EA_chr_file[i]);
			DNAStirrer ds = new DNAStirrer(afr, 1, 10000, AdmixtureConstant.Without_Genetic_Drift, w);
			ds.DNAStir(1);
			EA_DNAPool.add(ds);
			ChromosomeGenerator cg = new ChromosomeGenerator(ds.AncestrySNPFreqencyPanel(), ds.CurrSNPOrigine());
			cg.setSeed(seed + i);
			EA_CG.add(cg);
		}

		double[] EA_disease_rate = { 0.5 };
		GenerateColony EA_GC = new GenerateColony(N_phe, seed, disease_chr, EA_disease_rate,
				hs, EA_DNAPool, EA_CG, pg, AdmixtureConstant.free_recombination);

		// specific components
		// family
		int EA_N_Fam = 100;
		int EA_N_Kid = 2;
		int EA_N_aff_Kid = 1;
		QualityControl EA_qc = new QualityControl(EA_N_aff_Kid, AdmixtureConstant.No_selection);
		EA_GC.GenerateNewFamHab(EA_N_Fam, EA_N_Kid, EA_qc);

		int EA_N_case = 50;
		int EA_N_control = 50;
		QualityControl EA_qc_c = new QualityControl(EA_N_case, EA_N_control, AdmixtureConstant.No_selection);
		EA_GC.GenerateCCHab(EA_N_case + EA_N_control, 1, EA_qc_c);

		Print2File p2f = new Print2File();
		p2f.addColony(AA_GC);
		p2f.addColony(EA_GC);
		try {
			p2f.printGenotype2file("ped.txt", "phe.txt", !AdmixtureConstant.printAllele);
			p2f.printUnrelatedIndividual("PCA_geno.txt", "PCA_phe.txt", !AdmixtureConstant.printAllele, true);
		} catch (IOException e) {
			e.printStackTrace(System.err);
		}
	}
}
