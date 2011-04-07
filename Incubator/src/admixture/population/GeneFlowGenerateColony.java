package admixture.population;

import java.util.ArrayList;
import java.util.Random;

import admixture.population.genome.GeneFlow;
import admixture.population.genome.chromosome.ChromosomeGenerator;
import admixture.population.genome.chromosome.FamilyGenome;
import admixture.population.phenotype.FamilyPhenotype;
import admixture.population.phenotype.QualityControl;

public class GeneFlowGenerateColony extends GenerateColony {

	protected ArrayList<GeneFlow> GF;
	protected double[] pop_prop;
	protected Random rnd = new Random(2010);

	protected static abstract class Init<T extends Init<T>> extends GenerateColony.Init<T> {
		private ArrayList<GeneFlow> GF;
		private double[] pop_prop;

		public T GeneFlow(ArrayList<GeneFlow> gf) {
			this.GF = gf;
			return self();
		}

		public T popProportion(double[] pp) {
			this.pop_prop = new double[pp.length];
			System.arraycopy(pp, 0, pop_prop, 0, pp.length);
			for (int i = 1; i < pop_prop.length; i++) {
				pop_prop[i] = pop_prop[i] + pop_prop[i - 1];
			}
			return self();
		}

		public GeneFlowGenerateColony build() {
			return new GeneFlowGenerateColony(this);
		}
	}

	public static class Builder extends Init<Builder> {
		@Override
		protected Builder self() {
			return this;
		}
	}

	public GeneFlowGenerateColony(Init<?> init) {
		super(init);
		this.GF = init.GF;
		this.pop_prop = init.pop_prop;
		// TODO Auto-generated constructor stub
	}

	protected void generateFamilies(Habitat hab, int N_Fam, int N_Kid, QualityControl qc) {
		for (GeneFlow gf : GF) {
			gf.randomization();
		}
		for (int i = 0; i < N_Fam; i++) {
			FamilyGenome fg = new FamilyGenome(CurrFam + i + 1, N_Kid);
			FamilyPhenotype fp;
			int r = 0;
			do {
				for (int j = 0; j < GF.size(); j++) {
					if (r > 0 && control_chr == j) {
						continue;
					}
					int chrID = j;
					GeneFlow geneflow = GF.get(j);
					ChromosomeGenerator cg = ChrGenerator.get(j);
					hs.rev(geneflow.NumberOfSNP());
					int[][] f_g;
					int[][] m_g;
					if (j == control_chr) {
						int[][] f = geneflow.getAnIndividualInPool(i * 2);
						int[][] m = geneflow.getAnIndividualInPool(i * 2 + 1);
						f_g = new int[f.length][f[0].length];
						m_g = new int[f.length][f[0].length];
						System.arraycopy(f[0], 0, f_g[0], 0, f[0].length);
						System.arraycopy(f[1], 0, f_g[1], 0, f[1].length);
						System.arraycopy(m[0], 0, m_g[0], 0, m[0].length);
						System.arraycopy(m[1], 0, m_g[1], 0, m[1].length);						
					} else {
						int idx = 0;
						float f = rnd.nextFloat();
						while (f > pop_prop[idx])
							idx++;

						f_g = geneflow.sampleAnFounder(idx);

						f = rnd.nextFloat();
						idx = 0;
						while (f > pop_prop[idx])
							idx++;

						m_g = geneflow.sampleAnFounder(idx);
					}
					if (r == 0) {
						fg.addFamilyChromosome(cg.generateFamilySingleChromosome(chrID, f_g, m_g, N_Kid, hs,
								control_chr != j));
					} else {
						fg.setFamilyChromosome(j,
								cg.generateFamilySingleChromosome(chrID, f_g, m_g, N_Kid, hs, control_chr != j));
					}
				}
				if (!isNullHypothesis) {
					fp = pheGenerator.getGeneratePhenotypeAdmixtureLogistic(fg, disease_rate);
				} else {
					fp = pheGenerator.getGeneratePhenotypeAncestry(fg, disease_rate);
				}
				r++;
			} while (!qc.Accept(fp));
			hab.AddFamilyGenome(fg);
			hab.AddFamilyPhenotype(fp);
		}
		CurrFam += N_Fam;
	}
}
