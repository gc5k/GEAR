import java.util.Arrays;

import jsc.distributions.DiscreteUniform;

public class Test {

	public static void main(String[] args) {
		
		StringBuffer sb = new StringBuffer();
		for (int i = 0; i < 10; i++) {
			sb.append(i);
		}
		System.out.println(sb.toString());
		
		String[] S={"aa", "bb", "cc"};
		System.out.println(Arrays.binarySearch(S, "bb"));
	}
}
