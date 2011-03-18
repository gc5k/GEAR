package admixture.population;

import java.util.ArrayList;
import admixture.population.genome.chromosome.FamilyGenome;
import admixture.population.phenotype.FamilyPhenotype;
import arsenal.NewIt;
/**
*
* @author Guo-Bo Chen, chenguobo@gmail.com
*/
public class Habitat {

	private ArrayList<FamilyGenome> FamGenome;
	private int NumChr;
	private ArrayList<FamilyPhenotype> FamPhenotype;
	private int NumPhe;
	
	public Habitat() {
		FamGenome = NewIt.newArrayList();
		FamPhenotype = NewIt.newArrayList();
	}

	public void AddFamilyGenome(FamilyGenome fg) {
		FamGenome.add(fg);
	}

	public void InsertFamilyGenome(int idx, FamilyGenome fg) {
		FamGenome.set(idx, fg);
	}

	public void AddFamilyPhenotype(FamilyPhenotype fp) {
		FamPhenotype.add(fp);
	}
	
	public void InsertFamilyPhenotype(int idx, FamilyPhenotype fp) {
		FamPhenotype.add(fp);
	}
	
	public ArrayList<FamilyGenome> getFamilyGenome() {
		return FamGenome;
	}

	public ArrayList<FamilyPhenotype> getFamilyPhenotype() {
		return FamPhenotype;
	}	
}
