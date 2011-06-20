package admixture.population.phenotype;

import java.util.Arrays;
import java.util.Random;

import admixture.population.genome.chromosome.FamilyGenome;

public class PhenotypeGenerator {

	private boolean DEBUG = false;
	private int[] linked_chr;
	private int[] linked_loci;
	private String[] function_score;
	private double c_dev;
	private double c_effect;
	private double mu;
	private double[] gene_effect;
	private long seed = 2011;
	private Random rnd;

	public PhenotypeGenerator(String[] fs, double[] g_e, int[] c, int[] l, double u, double d_eff, double d_cov) {
		function_score = fs;
		gene_effect = g_e;
		linked_chr = c;
		linked_loci = l;
		mu = u;
		c_effect = d_eff;
		c_dev = d_cov;
		rnd = new Random(seed);
	}

	public void setFunction_score() {
	}

	public void setLikedChr(int[] c) {
		linked_chr = c;
	}

	public void setLikedLoci(int[] loci) {
		linked_loci = loci;
	}

	public void setSeed(long s) {
		rnd.setSeed(s);
	}

	public FamilyPhenotype getGeneratePhenotypeAdmixtureLogistic(FamilyGenome fg, double[] disease_rate) {
		// covariate is not considered
		fg.AscertainGenomeAncestry();
		double[][] p_p = new double[2][disease_rate.length];
		int[] p_s = new int[2];
		for (int i = 0; i < p_p.length; i++) {

			String g = fg.ParentGenotype(i, linked_chr, linked_loci);
			double logit = mu + rnd.nextGaussian() * c_dev * c_effect;
			int idx = Arrays.binarySearch(function_score, g);
			if (idx >= 0)
				logit += gene_effect[idx];

			double mix_prob = Math.exp(logit) / (1 + Math.exp(logit));

			double[] a = fg.ParentAncestry(i);
			for (int j = 0; j < a.length; j++) {
				mix_prob += a[j] * disease_rate[j];
				p_p[i][j] = a[j];
			}

			if (DEBUG) {
				System.out.println("Parent " + i + " : disease rate " + mix_prob);
			}

			p_s[i] = rnd.nextFloat() < mix_prob ? 1 : 0;
		}

		double[][] o_p = new double[fg.getNumberOffspring()][disease_rate.length];
		int[] o_s = new int[fg.getNumberOffspring()];
		for (int i = 0; i < o_p.length; i++) {
			String g = fg.OffspringGenotype(i, linked_chr, linked_loci);
			double logit = mu + rnd.nextGaussian() * c_dev * c_effect;
			int idx = Arrays.binarySearch(function_score, g);
			if (idx >= 0)
				logit += gene_effect[idx];

			double mix_prob = Math.exp(logit) / (1 + Math.exp(logit));

			double[] a = fg.OffspringAncestry(i);
			for (int j = 0; j < a.length; j++) {
				mix_prob += a[j] * disease_rate[j];
				o_p[i][j] = a[j];
			}

			if (DEBUG) {
				System.out.println("Parent " + i + " : disease rate " + mix_prob);
			}

			o_s[i] = rnd.nextFloat() < mix_prob ? 1 : 0;
		}

		FamilyPhenotype fp = new FamilyPhenotype(fg.getFamilyID(), p_p, p_s, o_p, o_s);
		return fp;
	}

	public FamilyPhenotype getGeneratePhenotypeAncestry(FamilyGenome fg, double[] disease_rate) {
		fg.AscertainGenomeAncestry();
		double[][] p_p = new double[2][disease_rate.length];
		int[] p_s = new int[2];
		for (int i = 0; i < p_p.length; i++) {
			double[] a = fg.ParentAncestry(i);
			double mix_prob = 0;
			for (int j = 0; j < a.length; j++) {
				mix_prob += a[j] * disease_rate[j];
				p_p[i][j] = a[j];
			}
			if (DEBUG) {
				System.out.println("Parent " + i + " : disease rate " + mix_prob);
			}
			p_s[i] = rnd.nextFloat() < mix_prob ? 1 : 0;
		}

		double[][] o_p = new double[fg.getNumberOffspring()][disease_rate.length];
		int[] o_s = new int[fg.getNumberOffspring()];
		for (int i = 0; i < o_p.length; i++) {
			double[] a = fg.OffspringAncestry(i);
			double mix_prob = 0;
			for (int j = 0; j < a.length; j++) {
				mix_prob += a[j] * disease_rate[j];
				o_p[i][j] = a[j];
			}
			if (DEBUG) {
				System.out.println("Kid " + i + " : disease rate " + mix_prob);
			}
			o_s[i] = rnd.nextFloat() < mix_prob ? 1 : 0;
		}

		FamilyPhenotype fp = new FamilyPhenotype(fg.getFamilyID(), p_p, p_s, o_p, o_s);
		return fp;
	}

}
