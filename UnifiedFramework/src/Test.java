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
		
		byte genotype1 = 1;
		genotype += (byte)(genotype1 << 2);
		byte genotype2 = 2;
		byte genotype3 = 3;
		genotype += (genotype2 << 4);
		genotype += (genotype3 << 6);

		byte mask = 3;
		byte g = -1;
		for(int i = 0; i < 8; i += 2) {
			byte gg = (byte) ((genotype>>i) & mask);
			System.out.println(i+ " gg " +gg);
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
