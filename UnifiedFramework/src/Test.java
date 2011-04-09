import java.util.Collections;
import java.util.HashMap;

import util.NewIt;


public class Test {
	public static void main(String[] args) {
        HashMap<String, Integer> modelCount = NewIt.newHashMap();
        modelCount.put("ss", new Integer(0));
        modelCount.put("aa", new Integer(12));
        modelCount.put("sss", new Integer(3));
        Integer v = Collections.max(modelCount.values());
        System.out.println(v);
		byte genotype = 3;
		genotype = (byte)(genotype << 2);
		genotype += 2;
		genotype = (byte)(genotype << 2);
		genotype += 3;
		genotype = (byte)(genotype << 2);
		genotype += 3;

		int mask = 3;
		byte g = -1;
		for(int i = 0; i < 8; i++) {
			int gg = (g>>i) & mask;
			System.out.println(i+ " " +gg);
		}

		System.out.println((int)genotype);
		byte m1 = 3;
		byte m2 = 12;
		byte m3 = 48;
		byte m4 = -128;
		byte g1 = (byte) (genotype & mask);
		byte g2 = (byte) ((genotype>>2) & mask);
		byte g3 = (byte) ((genotype>>4) & mask);
		byte g4 = (byte) ((genotype>>6) & mask);
		System.out.println(g1 + " " + g2 + " " + g3 + " " + g4);
	}
}
