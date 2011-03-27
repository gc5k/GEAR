package admixture.population;

import java.util.ArrayList;

import admixture.population.genome.DNAStirrer;
import admixture.population.genome.GeneFlow;
import admixture.population.genome.HotSpot;
import admixture.population.genome.chromosome.ChromosomeGenerator;
import admixture.population.phenotype.PhenotypeGenerator;

public class GeneFlowGenerateColony extends GenerateColony {

	protected ArrayList<GeneFlow> GF;
	public GeneFlowGenerateColony(Builder builder, ArrayList<GeneFlow> gf) {
		super(builder);
		GF = gf; 
		// TODO Auto-generated constructor stub
	}

	public GeneFlowGenerateColony(int np, long s, int dc, double[] dr, HotSpot h, ArrayList<DNAStirrer> dp, ArrayList<ChromosomeGenerator> cg,
			PhenotypeGenerator p, boolean rf, boolean isNull) {
		super(np, s, dc, dr, h, dp, cg, p, rf, isNull);
		// TODO Auto-generated constructor stub
	}

}
