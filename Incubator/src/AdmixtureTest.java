
import java.util.ArrayList;
import admixture.AdmixtureConstant;
import admixture.ChromosomeGenerator;
import admixture.DNAStirrer;
import admixture.Habitat;
import admixture.HotSpot;
import admixture.chromosome.FamilyGenome;
import admixture.chromosome.FamilySingleChromosome;

/**
 *
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */

public class AdmixtureTest {

	public static void main(String[] args) {
		
		ArrayList<DNAStirrer> DNAPool = new ArrayList<DNAStirrer>();
		ArrayList<ChromosomeGenerator> CG = new ArrayList<ChromosomeGenerator>();
		for (int i = 0; i < 2; i++) {
			DNAStirrer ds = new DNAStirrer(2, 10000, 5* (i+1), AdmixtureConstant.With_Genetic_Drift);
			ds.DNAStir(1);
			DNAPool.add(ds);
			ChromosomeGenerator cg = new ChromosomeGenerator(ds.SNPPanel(), ds.PostSNPAncestralProb());
			CG.add(cg);
		}

		Habitat Hab = new Habitat();
		for (int i = 0; i < 10; i++) {
			FamilyGenome fg = new FamilyGenome(i);
			for (int j = 0; j < 2; j++) {
				int num_kid = 2;
				int chrID = j;
				DNAStirrer ds = DNAPool.get(j);
				ChromosomeGenerator cg = CG.get(j);
				HotSpot hs = new HotSpot(ds.NumberOfSNP());
				hs.GenerateRecombination(AdmixtureConstant.free_recombination);
				fg.add(cg.generateFamilySingleChromosome(chrID, num_kid, hs.FatherRecombinationFraction(), hs.FatherRecombinationFraction()));
				fg.set(j, cg.generateFamilySingleChromosome(chrID, num_kid, hs.FatherRecombinationFraction(), hs.FatherRecombinationFraction()));
			}
			Hab.add(fg);
		}
		System.out.println();
	}
}
