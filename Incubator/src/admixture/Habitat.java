package admixture;

import java.util.ArrayList;
import admixture.chromosome.FamilyGenome;
import admixture.phenotype.FamilyPhenotype;
import arsenal.NewIt;
/**
*
* @author Guo-Bo Chen, chenguobo@gmail.com
*/
public class Habitat {

	private ArrayList<FamilyGenome> FamGenome;
	private ArrayList<FamilyPhenotype> FamPhenotype;

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
}
