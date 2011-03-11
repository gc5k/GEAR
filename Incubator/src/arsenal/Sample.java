package arsenal;

import java.util.Random;
import jsc.distributions.DiscreteUniform;
import org.apache.commons.lang3.ArrayUtils;
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

	public static int[] SampleIndex(int start, int end, int n, boolean replacement) {
		if(start > end) {
			return null;
		}
		int[] idx = new int[n];
		if(replacement) {
			for(int i = 0; i < idx.length; i++) {
				idx[i] = (int) U.nextInt(end - start + 1) + start;
			}
		} else if ( n == (end - start + 1)) {
			for(int i = 0; i < n; i++) {
				idx[i] = start + i;
			}
		} else {
			int[] url = new int[end - start + 1];
			for(int i = 0; i < url.length; i++) url[i] = i + start;
			for(int i = 0; i < n; i++) {
				int j = (int) U.nextInt(url.length);
				idx[i] = url[j];
				url = ArrayUtils.remove(url, j);
			}
		}
		return idx;
	}
}
