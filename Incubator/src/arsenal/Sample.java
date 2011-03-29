package arsenal;

import java.util.Random;

/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class Sample {
	public static long seed = 2011;
	public static Random U = new Random(2011);

	public static void setSeed(long s) {
		U.setSeed(s);
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
