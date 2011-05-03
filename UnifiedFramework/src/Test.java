import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.SortedMap;
import java.util.Random;
import java.util.TreeMap;
import java.util.Map.Entry;

import util.NewIt;

public class Test {
	private final static Comparator<byte[]> byteArrayComparator = new Comparator<byte[]>() {
		public int compare(final byte[] o1, final byte[] o2) {
			if (o1.length != o2.length) {
				throw new RuntimeException(
						"byte arrays are not the same length! first: "
								+ Arrays.toString(o1) + " second: " + Arrays.toString(o2));
			}
			int comparisonResult = 0;
			for (int index = 0; index < o1.length; ++index) {
				comparisonResult = o1[index] - o2[index];
				if (comparisonResult != 0) {
					break;
				}
			}
			return comparisonResult;
		}
	};
	public static void main(String[] args) {
		byte[][] d = new byte[1000][10];
		for(int i = 0; i < d.length; i++) {
			for(int j = 0; j < d[i].length; j++) {
				d[i][j] = (byte) j;
			}
		}
		long t = System.currentTimeMillis();
		
		SortedMap<byte[], Integer> m = new TreeMap<byte[], Integer>(byteArrayComparator);
		Random rand = new Random();
		for(int i = 0; i < d.length; i++) {
			final byte[] key = new byte[3];
			for(int j = 0; j < 3; j++) {
				key[j] = d[i][j];
			}
			Integer D = null;
			if(m.containsKey(key)) {
				D = m.get(key);
				D++;
			} else {
				D = new Integer(1);
			}
			m.put(key, D);
		}
		
		for(Entry<byte[], Integer> entry : m.entrySet()) {
			byte[] key = entry.getKey();
			Integer v = entry.getValue();
			System.out.println(key[0] + " " + key[1] + " " + key[2] + " " + v);
		}
	}
}
