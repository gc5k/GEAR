package admixture;

public class Tools {

	static public double[][] Generate_Matrix(double[] vector) {
		double[][] M = new double[vector.length][vector.length];
		for(int i = 0; i < vector.length; i++) {
			for(int j = 0; j < vector.length; j++) {
				M[i][j] = vector[i] * vector[j]; 
			}
		}
		return M;
	}
}
