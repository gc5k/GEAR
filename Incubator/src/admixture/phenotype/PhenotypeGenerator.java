package admixture.phenotype;

import java.util.Arrays;
import java.util.Random;

import admixture.chromosome.FamilyGenome;
import admixture.phenotype.FamilyPhenotype;

public class PhenotypeGenerator {
	public static int[] linked_chr;
	public static int[] linked_loci;
	public static String[] function_score;
	public static double dev_cov;
	public static double[] gene_effect;
	public static Random rnd = new Random();
	public static void setFunction_score(String[] fs, double[] g) {
		function_score = fs;
		gene_effect = g;
	}

	public static void setLikedChr(int[] c) {
		linked_chr = c;
	}

	public static void setLikedLoci(int[] loci) {
		linked_loci = loci;
	}

	public static FamilyPhenotype getGeneratePhenotype (FamilyGenome fg) {
		double[][] p_p = new double[2][1];
		int[] p_s = new int[2];
		for(int i = 0; i < p_p.length; i++) {
			String g = fg.ParentGenotype(i, linked_chr, linked_loci);
			int fun_idx = Arrays.binarySearch(function_score, g);
			p_p[i][0] = fun_idx>=0 ? rnd.nextGaussian() * dev_cov + gene_effect[fun_idx]:rnd.nextGaussian() * dev_cov;
			double prob = Math.exp(p_p[i][0]) / (1 + Math.exp(p_p[i][0]));
			p_s[i] = rnd.nextFloat()< prob ? 1:0;
		}

		double[][] o_p = new double[fg.getNumberOffspring()][1];
		int[] o_s = new int[fg.getNumberOffspring()];
		for(int i = 0; i < o_p.length; i++) {
			String g = fg.OffspringGenotype(i, linked_chr, linked_loci);
			int fun_idx = Arrays.binarySearch(function_score, g);
			p_p[i][0] = fun_idx>=0 ? rnd.nextGaussian() * dev_cov + gene_effect[fun_idx]:rnd.nextGaussian() * dev_cov;
			double prob = Math.exp(p_p[i][0]) / (1 + Math.exp(p_p[i][0]));
			p_s[i] = rnd.nextFloat()< prob ? 1:0;
		}
		FamilyPhenotype fp = new FamilyPhenotype(p_p, p_s, o_p, o_s);
		return fp;
	}
}
