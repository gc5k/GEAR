package gear.util;

import java.util.ArrayList;

//code is from https://github.com/hmkcode/Java/blob/master/java-combinations/Combination.java

public class Combination {

	/*
	 * public static void main(String[] args){ Integer[] elements = new
	 * Integer[] {1,2,3,4}; int[][] cL=combination(elements,3); }
	 */

	public static int[][] combination(Object[] elements, int K) 
	{
		ArrayList<int[]> comList = new ArrayList<int[]>();
		// get the length of the array
		// e.g. for {'A','B','C','D'} => N = 4
		int N = elements.length;

		if (K > N)
		{
			System.out.println("Invalid input, K > N");
			return null;
		}
		// calculate the possible combinations
		// e.g. c(4,2)
		c(N, K);

		// get the combination by index
		// e.g. 01 --> AB , 23 --> CD
		int combination[] = new int[K];

		// position of current index
		// if (r = 1) r*
		// index ==> 0 | 1 | 2
		// element ==> A | B | C
		int r = 0;
		int index = 0;

		while (r >= 0) 
		{
			// possible indexes for 1st position "r=0" are "0,1,2" --> "A,B,C"
			// possible indexes for 2nd position "r=1" are "1,2,3" --> "B,C,D"

			// for r = 0 ==> index < (4+ (0 - 2)) = 2
			if (index <= (N + (r - K))) 
			{
				combination[r] = index;

				// if we are at the last position print and increase the index
				if (r == K - 1)
				{

					// do something with the combination e.g. add to list or
					// print

//					print(combination, elements);
					int[] idx = new int[combination.length];
					System.arraycopy(combination, 0, idx, 0, combination.length);
					comList.add(idx);
					index++;
				}
				else
				{
					// select index for next position
					index = combination[r] + 1;
					r++;
				}
			}
			else 
			{
				r--;
				if (r > 0)
					index = combination[r] + 1;
				else
					index = combination[0] + 1;
			}
		}
		return comList.toArray(new int[0][]);
	}

	public static int c(int n, int r) 
	{

		int nf = fact(n);
		int rf = fact(r);
		int nrf = fact(n - r);
		int npr = nf / nrf;
		int ncr = npr / rf;

		return ncr;
	}

	public static int fact(int n) 
	{
		if (n == 0)
			return 1;
		else
			return n * fact(n - 1);
	}

	// public static void print(int[] combination, Object[] elements){
	//
	// String output = "";
	// for(int z = 0 ; z < combination.length;z++){
	// output += elements[combination[z]];
	// }
	// System.out.println(output);
	// }

}