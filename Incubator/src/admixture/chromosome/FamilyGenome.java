package admixture.chromosome;

import java.util.ArrayList;

public class FamilyGenome extends ArrayList<FamilySingleChromosome> {
	private static final long serialVersionUID = 1L;
	private int FamID;
//	private ArrayList<FamilySingleChromosome> FG = new ArrayList<FamilySingleChromosome>();

	public FamilyGenome(int fid) {
		FamID = fid;
	}
}
