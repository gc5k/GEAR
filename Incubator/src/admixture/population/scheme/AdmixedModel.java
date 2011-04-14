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
		Parameter p = new Parameter();
		p.commandListenor(args);

		long seed = p.seed;
		int control_chr = p.control_chr;
		boolean isNullHypothesis = p.isNullHypothesis;

		// specific components
		// family
		int N_Fam = p.family[0];
		int N_Kid = p.kid[0];
		int N_aff_Kid = p.affectedKid[0];

		int N_case = p.cases[0];
		int N_control = p.controls[0];

		// logistic regression
		String[] f = p.genotypeFunction;
		double[] g_e = p.genotypeEffect;
		int[] chr = p.diseaseChr;
		int[] loci = p.diseaseLocus;
		double mu = p.mu;
		double dev = p.covariateEffect;

		double[] prevalence = p.popPrevalence;
		int N_phe = 3;

		String[] chr_file = p.AIM_file;
		int[] AIM_number = p.aim;
		double[] pop_proportion = p.popProportion;

		PhenotypeGenerator pg = new PhenotypeGenerator(f, g_e, chr, loci, mu, dev);
		HotSpot hs = new HotSpot();

		ArrayList<DNAStirrer> DNAPool = new ArrayList<DNAStirrer>();
		ArrayList<GeneFlow> GF = new ArrayList<GeneFlow>();
		ArrayList<ChromosomeGenerator> CG = new ArrayList<ChromosomeGenerator>();
		for (int i = 0; i < chr_file.length; i++) {
			AlleleFrequencyReader afr = new AlleleFrequencyReader(chr_file[i], AIM_number[i]);
			DNAStirrer ds = new DNAStirrer(afr, 1, 10000, AdmixtureConstant.Without_Genetic_Drift, pop_proportion);
			ds.DNAStir(1);
			DNAPool.add(ds);
			GeneFlow gf = new GeneFlow(afr, 1000, pop_proportion, hs);
			gf.mating(p.generation);
			GF.add(gf);
			ChromosomeGenerator cg = new ChromosomeGenerator(afr.getAlleleFreq(), pop_proportion);
			cg.setSeed(seed + i);
			CG.add(cg);
		}

		GeneFlowGenerateColony GC = new GeneFlowGenerateColony.Builder().DNAPool(DNAPool).diesaseRate(prevalence)
				.hotSpot(hs).popProportion(pop_proportion).numPhenotype(N_phe).ChrGenerator(CG).PheGenerator(pg).GeneFlow(GF).seed(
						seed).diseaseChr(control_chr).isNullHypothesis(isNullHypothesis).build();

		for (int rep = 0; rep < 2; rep++) {
			QualityControl qc = new QualityControl(N_aff_Kid, AdmixtureConstant.FamilyExactAffected);
			QualityControl qc_c = new QualityControl(N_case, N_control, AdmixtureConstant.CaseControl);
			GeneFlowGenerateColony.setCurrFamilyID(0);
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
