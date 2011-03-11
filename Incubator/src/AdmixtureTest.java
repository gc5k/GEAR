
import java.util.ArrayList;
import admixture.AdmixtureConstant;
import admixture.ChromosomeGenerator;
import admixture.DNAStirrer;
import admixture.Habitat;
import admixture.HotSpot;
import admixture.chromosome.FamilyGenome;
import admixture.phenotype.FamilyPhenotype;

import admixture.phenotype.PhenotypeGenerator;
/**
 *
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */

public class AdmixtureTest {

	public static void main(String[] args) {

		int disease_chr = -1;
		String[] f = {"0000", "0101", "1111"};
		double[] g_e = {1, 0.5, 1};
		int[] chr = {1, 1};
		int[] loci = {0, 1};
		double mu =  0;
		double dev = Math.sqrt(10);
		PhenotypeGenerator pg = new PhenotypeGenerator(f, g_e, chr, loci, mu, dev);

		Habitat Hab = new Habitat();
		double[] disease_pop = {0.2, 0.5};

		ArrayList<DNAStirrer> DNAPool = new ArrayList<DNAStirrer>();
		ArrayList<ChromosomeGenerator> CG = new ArrayList<ChromosomeGenerator>();
		for (int i = 0; i < 2; i++) {
			DNAStirrer ds = new DNAStirrer(2, 10000, 5 * (i+1), AdmixtureConstant.Without_Genetic_Drift);
			ds.DNAStir(1);
			DNAPool.add(ds);
			ChromosomeGenerator cg = new ChromosomeGenerator(ds.SNPPanel());
			CG.add(cg);
		}
		HotSpot hs = new HotSpot(2010);
		for (int i = 0; i < 1; i++) {
			int num_kid = 2;
			FamilyGenome fg = new FamilyGenome(i, num_kid);

			for (int j = 0; j < 2; j++) {
				int chrID = j;
				DNAStirrer ds = DNAPool.get(j);
				ChromosomeGenerator cg = CG.get(j);
				hs.rev(ds.NumberOfSNP());
				hs.GenerateRecombination(AdmixtureConstant.free_recombination);
				int[] f_hotspot = hs.getHotSpot();
				hs.GenerateRecombination(AdmixtureConstant.free_recombination);
				int[] m_hotspot = hs.getHotSpot();

				fg.addFamilyChromosome(cg.generateFamilySingleChromosome(chrID, num_kid, f_hotspot, m_hotspot, ds.PostSNPAncestralProb(), disease_chr ==j));
			}
			fg.AscertainGenomeAncestry();
//			fg.printAncestry();
//			fg.printGenome();
			FamilyPhenotype fp = pg.getGeneratePhenotypeAncestry(fg, disease_pop);
//			fp.print();
			Hab.AddFamilyGenome(fg);
			Hab.AddFamilyPhenotype(fp);
		}
		System.out.println();
	}
}
