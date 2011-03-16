package admixture;

import java.util.Random;

import admixture.chromosome.FamilySingleChromosome;

/**
*
* @author Guo-Bo Chen, chenguobo@gmail.com
*/

public class ChromosomeGenerator {
	private boolean DEBUG = false;
	private int N_snp;
	private Random rnd;
	private double[] snp_panel;
	private double[] LD; // reserve for the future
	private int[][] hotspot; // there two kinds of recombination fractions
	// 1: free of recombination that each element equals 0.5, hotspot=[0.5,0.5,0.5,...]
	// it is generated in the method recombinationFree()
	// 2: cross over between snps: hotspot looks like=[0, 0, 0, 1, 0, 0, 1, ...]
	// it is generated in the method recombination()
	
	// Note: hotspot[0] for paternal recombination, hotspot[1] for maternal recombination;

	public ChromosomeGenerator(double[] sp) {
		snp_panel = sp;
		rnd = new Random(2011);
		hotspot = new int[2][];
		N_snp = sp.length;
	}

	public FamilySingleChromosome generateFamilySingleChromosome(int cID, int k, int[] father_hotspot, int[] mother_hotspot, double[][][] post_snp_ancestry, boolean disease_linked) {
		hotspot[0] = father_hotspot;
		hotspot[1] = mother_hotspot;		
		int[][][] pg = new int[2][][];
		pg[0] = generateFounderChr(AdmixtureConstant.Without_LD);
		pg[1] = generateFounderChr(AdmixtureConstant.Without_LD);

		int[][][] og = new int[k][][];
		for (int i = 0; i < k; i++) {
			og[i] = generateOffspringChr(pg);
		}
		FamilySingleChromosome fsc = new FamilySingleChromosome(cID, pg, og, disease_linked);
		if(post_snp_ancestry != null) {
			fsc.AscertainParentSingleChromosomeAncestry(post_snp_ancestry);
			fsc.AscertainOffspringSingleChromosomeAncestry(post_snp_ancestry);
		}
		return fsc;
	}

	private int[][] generateFounderChr(boolean ld) {
		int[][] diploid = new int[2][N_snp];
		for (int i = 0; i < 2; i++) {
			for (int k = 0; k < N_snp; k++) {
				if(!ld) {
					diploid[i][k] = rnd.nextFloat() < snp_panel[k] ? 0:1;
				} else {
					// when there is LD pattern;
				}
			}
		}
		return diploid;
	}

	private int[][] generateOffspringChr(int[][][] pg) {//pg[0] for paternal chromosome, pg[1] for maternal.
		int[][] diploid = new int[2][N_snp];
		for (int i = 0; i < hotspot.length; ++i) {
			int chromatid = rnd.nextBoolean() ? 0 : 1;
//			System.out.println("-----"+ chromatid);
//			for (int j = 0; j < hotspot[i].length; j++) {
//				System.out.println(" =" + hotspot[i][j]);
//			}
			for (int j = 0; j < hotspot[i].length - 1; j++) {
				chromatid = 1 - chromatid;
				for (int k = hotspot[i][j]; k <= hotspot[i][j+1]; k++) {
					diploid[i][k] = pg[i][chromatid][k];
				}
			}
		}
		return diploid;
	}

	public void setSeed(long s) {
		rnd.setSeed(s);
	}

	public static void main(String[] args) {
		double[] snp_freq = { 0.5, 0.5, 0.1, 0.1, 0.1 };
		HotSpot hs = new HotSpot();
		hs.rev(snp_freq.length);
		hs.GenerateRecombination(AdmixtureConstant.free_recombination);
	}
}
