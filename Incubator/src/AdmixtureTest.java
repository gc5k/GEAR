
import java.util.ArrayList;
import admixture.AdmixtureConstant;
import admixture.ChromosomeGenerator;
import admixture.DNAStirrer;
import admixture.Habitat;
import admixture.HotSpot;
import admixture.chromosome.FamilyGenome;

/**
 *
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */

public class AdmixtureTest {

	public static void main(String[] args) {

		Habitat Hab = new Habitat();

		ArrayList<DNAStirrer> DNAPool = new ArrayList<DNAStirrer>();
		ArrayList<ChromosomeGenerator> CG = new ArrayList<ChromosomeGenerator>();
		for (int i = 0; i < 2; i++) {
			DNAStirrer ds = new DNAStirrer(2, 10000, 10 * (i+1), AdmixtureConstant.With_Genetic_Drift);
			ds.DNAStir(1);
			DNAPool.add(ds);
			ChromosomeGenerator cg = new ChromosomeGenerator(ds.SNPPanel());
			CG.add(cg);
		}

		for (int i = 0; i < 10; i++) {
			int num_kid = 2;
			FamilyGenome fg = new FamilyGenome(i, num_kid);

			for (int j = 0; j < 2; j++) {
				int chrID = j;
				DNAStirrer ds = DNAPool.get(j);
				ChromosomeGenerator cg = CG.get(j);

				HotSpot hs = new HotSpot(ds.NumberOfSNP());
				hs.GenerateRecombination(AdmixtureConstant.free_recombination);
				int[] f_hotspot = hs.getHotSpot();
				hs.GenerateRecombination(AdmixtureConstant.free_recombination);
				int[] m_hotspot = hs.getHotSpot();

				fg.add(cg.generateFamilySingleChromosome(chrID, num_kid, f_hotspot, m_hotspot, ds.PostSNPAncestralProb()));
			}
			fg.AscertainGenomeAncestry();
			fg.printAncestry();
			Hab.AddFamilyGenome(fg);
		}
		System.out.println();
	}
}
