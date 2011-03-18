package admixture.population.genome;

import java.util.Arrays;

import admixture.AdmixtureConstant;
import arsenal.Sample;

import jsc.distributions.Poisson;
import jsc.distributions.Binomial;

/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */

public class HotSpot {

	private boolean DEBUG = false;
	private double len_Morgan = 1; // length of the chromosome measured by Morgan
	private int N_snp = 1;
	private long seed = 2010;
	private Binomial B;
	private Poisson Po;

	private int[] rec_hotspot;

	public HotSpot() {
		B = new Binomial(N_snp, 0.5);
		Po = new Poisson(len_Morgan);	
		B.setSeed(seed);
		Po.setSeed(seed);		
	}

	public void rev(int N_s) {
		N_snp = N_s;
		len_Morgan = 1;
		
		B.setN(N_s);	
	}
	
	public void GenerateRecombination(boolean rf) {
		if (rf) {
			recombinationFree();
		} else {
			recombination();
		}
	}

	public long getSeed() {
		return seed;
	}

	public void setSeed(long s) {
		seed = s;
		B.setSeed(s);
		Po.setSeed(s);
	}

	private void recombinationFree() {

		int crossover;
		do {
			crossover = (int) B.random() - 1;
		} while (crossover < 0);
//		System.out.println("crossover :" + crossover);
		rec_hotspot = new int[crossover + 2];
		if (crossover > 0) {
			int[] hs = Sample.SampleIndex(1, N_snp - 1, crossover, AdmixtureConstant.Without_replacement);
			System.arraycopy(hs, 0, rec_hotspot, 1, hs.length);
		}
		rec_hotspot[0] = 0;
		rec_hotspot[rec_hotspot.length - 1] = N_snp - 1;
		Arrays.sort(rec_hotspot);

		if (DEBUG) print();
	}

	private void recombination() {

		Poisson p = new Poisson(len_Morgan);

		int crossover = 0;
		do {
			crossover = (int) p.random();
		} while (crossover >= N_snp - 2);

		rec_hotspot = new int[crossover + 2];
		if (crossover > 0) {
			int[] hs = Sample.SampleIndex(1, N_snp, crossover, AdmixtureConstant.Without_replacement);
			System.arraycopy(hs, 0, rec_hotspot, 1, hs.length);
		}
		Arrays.sort(rec_hotspot);

		if (DEBUG)
			print();
	}

	public int[] getHotSpot() {
		return rec_hotspot;
	}


	private void print() {
		if (DEBUG) {
			for (int i = 0; i < rec_hotspot.length; i++) {
				System.out.print(rec_hotspot[i] + " ");
			}
			System.out.println();
		}
	}

	public static void main(String[] args) {
		double[] snp_freq = { 0.5, 0.5, 0.1, 0.1, 0.1 };
		HotSpot hs = new HotSpot();
		hs.rev(snp_freq.length);
		hs.GenerateRecombination(AdmixtureConstant.free_recombination);
	}
}
