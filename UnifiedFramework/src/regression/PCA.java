package regression;

import java.util.Random;

import javastat.util.BasicStatistics;
import javastat.util.DataManager;
import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import util.DataOperator;

/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class PCA {
	private double[][] principalComponents;
	private double[] variance;
	private double[][] data;
	private double[][] covariance;

	public PCA(double[][] d) {
		data = d;
	}

	public double[][] principalComponent() {
		BasicStatistics basicStatistics = new BasicStatistics();
		covariance = basicStatistics.correlationMatrix(data);
		System.out.println("cov");
		for(int i = 0; i < covariance.length; i++) {
			for (int j = 0; j < covariance[i].length; j++) {
				System.out.print(covariance[i][j] + "\t");
			}
			System.out.println();
		}
		
		Matrix covarianceMatrix = new Matrix(covariance);
		EigenvalueDecomposition eigenCovMatrix = new EigenvalueDecomposition(
				covarianceMatrix);
		double[] eigenvalueCovMatrix = eigenCovMatrix.getRealEigenvalues();

		double[][] eigenvectorCov = eigenCovMatrix.getV().getArray();
		DataManager md = new DataManager();
		md.dataSort(eigenvalueCovMatrix, eigenvectorCov);
		System.out.println("eigenvalue");
		for (int i = 0; i < eigenvalueCovMatrix.length; i++) {
			System.out.println(eigenvalueCovMatrix[i]);
		}
		principalComponents = new double[data.length][data.length];
		variance = new double[data.length];
		for (int i = data.length - 1; i >= 0; i--) {
			variance[data.length - 1 - i] = eigenvalueCovMatrix[i];
			for (int k = 0; k < data.length; k++) {
				principalComponents[data.length - 1 - i][k] = eigenvectorCov[k][i];
			}
		}
		return principalComponents;
	}

	public double[] PCscore(int idx) {
		double[] ps = new double[data[idx].length];
		double[][] cd = DataOperator.Normalize(data, true, true);
		for(int i = 0; i < cd[idx].length; i++) {
			for(int j = 0; j < cd.length; j++) {
				ps[i] += cd[j][i]*principalComponents[idx][j];
			}
		}
		return ps;
	}

	public static void main (String[] args) {
		Random rnd = new Random(10);
		double[][] data = new double[2][50];
		for(int i = 0; i < data.length; i++) {
			for(int j = 0; j < data[i].length; j++) {
				if (i == 0) {
					data[i][j] = rnd.nextBoolean() ? 1:0;
				} else {
					data[i][j] = data[i-1][j]; 
				}
				System.out.print(data[i][j] + "\t");
			}
			System.out.println();
		}
		System.out.println();
		double[][] nd = DataOperator.Normalize(data, true, true);
		for(int i = 0; i < nd.length; i++) {
			for(int j = 0; j < nd[i].length; j++) {
				data[i][j] = rnd.nextBoolean() ? 1:0;
				System.out.print(nd[i][j] + "\t");
			}
			System.out.println();
		}

		PCA pca = new PCA(data);

		double[][] pc = pca.principalComponent();
		System.out.println("eigen vectors");
		for(int i = 0; i < pc.length; i++) {
			for(int j = 0; j < pc[i].length; j++) {
				System.out.print(pc[i][j] + "\t");
			}
			System.out.println();
		}
		System.out.println("PC score");
		double[] ps = pca.PCscore(0);
		for(int i = 0; i < ps.length; i++) {
			System.out.println(ps[i]);
		}
	}
}
