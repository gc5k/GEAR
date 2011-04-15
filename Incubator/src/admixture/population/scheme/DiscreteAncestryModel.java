package admixture.population.scheme;

import java.io.File;
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

		Parameter p = new Parameter();
		p.commandListenor(args);
		String unified = "unified";
		File file = new File(p.dir + unified);
		file.mkdir();
		String dir_unified = p.dir + unified + System.getProperty("file.separator");
		
		String fam = "family";
		file = new File(p.dir + fam);
		file.mkdir();
		String dir_fam = p.dir + fam + System.getProperty("file.separator");
		
		String c_c = "cc";
		file = new File(p.dir + c_c);
		file.mkdir();
		String dir_cc = p.dir + c_c + System.getProperty("file.separator");
		
		long seed = p.seed;
		int control_chr = p.control_chr;
		boolean isNullHypothesis = p.isNullHypothesis;
		// logistic regression

		String[] f = p.genotypeFunction;
		double[] g_e = p.genotypeEffect;
		int[] chr = p.diseaseChr;
		int[] loci = p.diseaseLocus;
		double mu = p.mu;
		double dev = p.covariateEffect;
		PhenotypeGenerator pg = new PhenotypeGenerator(f, g_e, chr, loci, mu, dev);
		HotSpot hs = new HotSpot();

		int N_phe = 2;

		// specify AA parameters
		double[] AA_disease_rate = new double[]{p.popPrevalence[0]};
		int AA_N_Fam = p.family[0];
		int AA_N_Kid = p.kid[0];
		int AA_N_aff_Kid = p.affectedKid[0];

		int AA_N_case = p.cases[0];
		int AA_N_control = p.controls[0];

		String[] AA_chr_file = p.AIM_file[0];
		int[] AIM_number = p.aim;
		double[] AA_w = new double[] { p.popProportion[0] };
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
		double[] EA_disease_rate = new double[]{ p.popPrevalence[1] };
		int EA_N_Fam = p.family[1];
		int EA_N_Kid = p.kid[1];
		int EA_N_aff_Kid = p.affectedKid[1];

		int EA_N_case = p.cases[1];
		int EA_N_control = p.controls[1];

		// generate EA
		String[] EA_chr_file = p.AIM_file[1];
		double[] w = new double[]{p.popProportion[1]};
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
				StringBuilder Gsb = new StringBuilder(rep + "ped.txt");
				Gsb.insert(0, dir_unified);
				StringBuilder Psb = new StringBuilder(rep + "phe.txt");
				Psb.insert(0, dir_unified);
				StringBuilder L_Gsb = new StringBuilder(rep + "L_ped.txt");
				L_Gsb.insert(0, dir_unified);
				StringBuilder L_Psb = new StringBuilder(rep + "L_phe.txt");
				L_Psb.insert(0, dir_unified);
				p2f.printGenotype2file(Gsb.toString(), Psb.toString(), !AdmixtureConstant.printAllele,
						!printLinked);
				p2f.printGenotype2file(L_Gsb.toString(), L_Psb.toString(), AdmixtureConstant.printAllele,
						printLinked);

				StringBuilder F_Gsb = new StringBuilder(rep + "ped.txt");
				F_Gsb.insert(0, dir_fam);
				StringBuilder F_Psb = new StringBuilder(rep + "phe.txt");
				F_Psb.insert(0, dir_fam);
				StringBuilder L_F_Gsb = new StringBuilder(rep + "L_ped.txt");
				L_F_Gsb.insert(0, dir_fam);
				StringBuilder L_F_Psb = new StringBuilder(rep + "L_phe.txt");
				L_F_Psb.insert(0, dir_fam);
				p2f.printFamilyGenotype2file(F_Gsb.toString(), F_Psb.toString(), !AdmixtureConstant.printAllele,
						!printLinked);
				p2f.printFamilyGenotype2file(L_F_Gsb.toString(), L_F_Psb.toString(), AdmixtureConstant.printAllele,
						printLinked);

				StringBuilder C_Gsb = new StringBuilder(rep + "ped.txt");
				C_Gsb.insert(0, dir_cc);
				StringBuilder C_Psb = new StringBuilder(rep + "phe.txt");
				C_Psb.insert(0, dir_cc);
				StringBuilder L_C_Gsb = new StringBuilder(rep + "L_ped.txt");
				L_C_Gsb.insert(0, dir_cc);
				StringBuilder L_C_Psb = new StringBuilder(rep + "L_phe.txt");
				L_C_Psb.insert(0, dir_cc);
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
