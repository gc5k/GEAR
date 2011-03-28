package admixture.population.scheme;

import java.io.IOException;
import java.util.ArrayList;

import admixture.AdmixtureConstant;
import admixture.population.AlleleFrequencyReader;
import admixture.population.GeneFlowGenerateColony;
import admixture.population.genome.DNAStirrer;
import admixture.population.genome.GeneFlow;
import admixture.population.genome.HotSpot;
import admixture.population.genome.chromosome.ChromosomeGenerator;
import admixture.population.phenotype.PhenotypeGenerator;
import admixture.population.phenotype.QualityControl;

public class AdmixedModel {
	public static void main(String[] args) {

		long seed = 2012;
		int control_chr = 0;
		boolean isNullHypothesis = false;
		
		// specific components
		// family
		int N_Fam = 300;
		int N_Kid = 2;
		int N_aff_Kid = 1;

		int N_case = 100;
		int N_control = 100;

		// logistic regression
		String[] f = { "0000", "1010", "1111" };
		double[] g_e = { 1, 0.5, 1 };
		int[] chr = { 1, 1 };
		int[] loci = { 0, 1 };
		double mu = 0;
		double dev = Math.sqrt(10);

		double[] disease_rate = { 0, 0 };
		int N_phe = 3;

		String[] chr_file = {"allele_freq_chr1.txt", "allele_freq_chr2.txt"};
		double[] w = {0.8, 0.2};

		PhenotypeGenerator pg = new PhenotypeGenerator(f, g_e, chr, loci, mu, dev);
		HotSpot hs = new HotSpot();
		
		ArrayList<DNAStirrer> DNAPool = new ArrayList<DNAStirrer>();
		ArrayList<GeneFlow> GF = new ArrayList<GeneFlow>();
		ArrayList<ChromosomeGenerator> CG = new ArrayList<ChromosomeGenerator>();
		for (int i = 0; i < chr_file.length; i++) {
			AlleleFrequencyReader afr = new AlleleFrequencyReader(chr_file[i]);
			DNAStirrer ds = new DNAStirrer(afr, 1, 10000, AdmixtureConstant.Without_Genetic_Drift, w);
			ds.DNAStir(1);
			DNAPool.add(ds);
			GeneFlow gf = new GeneFlow(afr, 1000, w, hs);
			if(i==0) gf.mating(2);
			GF.add(gf);
			ChromosomeGenerator cg = new ChromosomeGenerator(afr.getAlleleFreq(), w);
			cg.setSeed(seed + i);
			CG.add(cg);
		}

		GeneFlowGenerateColony GC = new GeneFlowGenerateColony.Builder().DNAPool(DNAPool).diesaseRate(disease_rate).hotSpot(hs).popProportion(w)
										.numPhenotype(N_phe).ChrGenerator(CG).PheGenerator(pg).GeneFlow(GF)
								.seed(seed).diseaseChr(control_chr).isNullHypothesis(isNullHypothesis).build();

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
		} catch (IOException e) {
			e.printStackTrace(System.err);
		}
	}
}
