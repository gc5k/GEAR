package arsenal;

/**
*
* @author Guo-Bo Chen, chenguobo@gmail.com
*/
public class Tools {

	public static double[][] Generate_Matrix(double[] vector) {
		double[][] M = new double[vector.length][vector.length];
		for(int i = 0; i < vector.length; i++) {
			for(int j = 0; j < vector.length; j++) {
				M[i][j] = vector[i] * vector[j]; 
			}
		}
		return M;
	}
	
	public static void Fill_Matrix(double[] scr, int[] idx, double v) {

		for(int i = 0; i < idx.length; i++) {
			scr[idx[i]] = v;
		}
	}

}
