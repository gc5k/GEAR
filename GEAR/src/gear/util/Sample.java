package gear.util;

import java.util.Random;

import parameter.Parameter;

/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */

public class Sample {
	public static long seed = 2011;
	public static Random U = new Random(Parameter.INSTANCE.seed);

	public static void setSeed(long s) {
		U.setSeed(s);
	}

	public static int[] sample(int[] s) {
		int[] seq = new int[s.length];
		System.arraycopy(s, 0, seq, 0, s.length);
		int[] index = new int[s.length];
		int len = s.length;
		for(int i = 0; i < index.length; i++) {
			int j = U.nextInt(len);
			len--;
			index[i] = seq[j];
			seq[j] = seq[len];
			seq[len] = index[i];
		}
		return index;
	}

	public static int[] SampleIndex(int start, int end, int n) {
		int[] seq = new int[end - start + 1];
		int[] index = new int[n];
		for (int i = 0; i < seq.length; i++)
			seq[i] = i + start;
		int len = seq.length;
		for(int i = 0; i < index.length; i++) {
			int j = U.nextInt(len);
			len--;			
			index[i] = seq[j];
			seq[j] = seq[len];
			seq[len] = index[i];
		}
		return index;
	}
	
	public static int[] SampleIndexWithReplacement(int start, int end, int n) {
		int[] seq = new int[end - start + 1];
		int[] index = new int[n];
		for (int i = 0; i < seq.length; i++)
			seq[i] = i + start;
		int len = seq.length;
		for(int i = 0; i < index.length; i++) {
			int j = U.nextInt(len);
			index[i] = seq[j];
		}
		return index;
	}
	

	public static void main(String[] args) {
		for (int i = 0; i < 10000; i++) {
			System.out.println(i);
			SampleIndex(0, 10000, 50);
		}
	}
}
