package score;

import java.util.*;
import org.apache.commons.math.linear.*;

public class MultiClassLogisticRegression {
	int category;
	double[][] Y;
	double[][] P;
	double[][] X;
	double[] B;
	double threshold = 0.001;
	double LogLikelihood_old;
	double LogLikelihood_new;

	public MultiClassLogisticRegression(int[] y, double[][] x, int k) {
		category = k;
		Y = new double[y.length][k - 1];
		P = new double[y.length][k - 1];
		for (int i = 0; i < y.length; i++) {
			if (y[i] < (k - 1)) {
				Y[i][y[i]] = 1;
			}
		}
		X = new double[x.length * (k - 1)][x[0].length * (k - 1)];
		for (int i = 0; i < k - 1; i++) {
			for (int j = 0; j < x.length; j++) {
				System.arraycopy(x[j], 0, X[i * x.length + j], x[j].length * i, x[j].length);
			}
		}
		B = new double[x[0].length * (k - 1)];
	}

	public void MLE() {
		RealMatrix Matrix_X = new RealMatrixImpl(X);
		RealMatrix Matrix_XT = Matrix_X.transpose();
		int iteration = 20;
		int iter = 0;
		RealMatrix B_old = new RealMatrixImpl(B);
		RealMatrix B_new = new RealMatrixImpl(B);
		LogLikelihood_old = convergence(B_old);
		do {
			B_old = B_new;
			LogLikelihood_old = LogLikelihood_new;
			B = B_old.getColumn(0);
			RealMatrix Matrix_W = new RealMatrixImpl(getWMatrix(B));
			RealMatrix Matrix_XT_W = Matrix_XT.multiply(Matrix_W);
			RealMatrix Matrix_XT_W_X = Matrix_XT_W.multiply(Matrix_X);
			RealMatrix Inv_XT_W_X = Matrix_XT_W_X.inverse();
			RealMatrix Inv_XT_W_X_XT = Inv_XT_W_X.multiply(Matrix_XT);
			RealMatrix Vector_Pre = new RealMatrixImpl(predicted(B));
			RealMatrix Vector_H = Inv_XT_W_X_XT.multiply(Vector_Pre);
			B_new = B_old.add(Vector_H);
			LogLikelihood_new = convergence(B_new);
			System.out.println(iter + "--->" + B_new);
			System.out.println(LogLikelihood_old - LogLikelihood_new);
		} while (iter++ < iteration && Math.abs(LogLikelihood_old - LogLikelihood_new) > threshold);
		System.out.println();
	}

	public void calculate_P(double[] b) {
		for (int i = 0; i < P.length; i++) {
			double[] p_ele = new double[category - 1];
			double sum = 0;
			for (int j = 0; j < P[i].length; j++) {
				double s = 0;
				for (int k = P.length * j; k < P.length * j + P.length; k++) {
					s += X[P.length * j][k] * b[k];
				}
				p_ele[j] = Math.exp(s);
				sum += p_ele[j];
			}
			for (int j = 0; j < P[i].length; j++) {
				P[i][j] = p_ele[j] / (1 + sum);
			}
		}
	}

	public double[][] predicted(double[] b) {
		double[][] res = new double[Y.length][Y[0].length];
		for (int i = 0; i < res.length; i++) {
			for (int j = 0; j < res[i].length; j++) {
				res[i][j] = Y[i][j] - P[i][j];
			}
		}
		return res;
	}

	public double[][] getWMatrix( double[] b) {
		double[][] W = new double[Y.length * (category - 1)][Y.length * (category - 1)];
		for (int i = 0; i < category - 1; i++) {
			for (int j = 0; j < category - 1; j++) {
				for (int k = 0; k < P.length; k++) {
					W[i * P.length + k][j * P.length + k] = i == j ? P[k][i] * (1 - P[k][j]) : P[k][i] * P[k][j] * (-1);
				}
			}
		}
		return W;
	}

	public boolean stop(RealMatrix b_old, RealMatrix b_new) {
		boolean stopit = true;
		for (int i = 0; i < b_old.getRowDimension(); i++) {
			if (Math.abs(b_old.getEntry(i, 0) - b_new.getEntry(i, 0)) > threshold) {
				stopit = false;
				break;
			}
		}
		return stopit;
	}

	public double convergence(RealMatrix b) {
		boolean converg = true;
		double L = 0;
		for (int i = 0; i < P.length; i++) {
			double sum_y = 0;
			double sum_p = 0;
			for (int j = 0; j < P[i].length; j++) {
				L += Y[i][j] * Math.log10(P[i][j]);
				sum_y += Y[i][j];
				sum_p += P[i][j];
			}
			L += (1-sum_y) * Math.log10(1-sum_p);
		}
		return L * (-1);
	}

	public static void main(String[] args) {
		int[] Y = { 2, 1, 2, 1, 1, 0, 2, 0, 0, 0 };
		double[][] X = { { 1, 15, 4 }, { 1, 30, 14 }, { 1, 31, 16 }, { 1, 31, 11 }, { 1, 32, 17 }, { 1, 29, 10 },
				{ 1, 30, 8 }, { 1, 31, 12 }, { 1, 32, 6 }, { 1, 40, 7 } };
		MultiClassLogisticRegression MCLR = new MultiClassLogisticRegression(Y, X, 3);
		MCLR.MLE();
	}
}
