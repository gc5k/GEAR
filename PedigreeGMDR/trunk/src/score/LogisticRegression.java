package score;

import java.util.*;
import org.apache.commons.math.linear.*;

public class LogisticRegression {
	double[] Y;
	double[] P;
	double[][] X;
	double[] B;
	double threshold = 0.001;
	double LogLikelihood_old;
	double LogLikelihood_new;

	LogisticRegression(double[] y, double[][] x) {
		Y = new double[y.length];
		System.arraycopy(y, 0, Y, 0, Y.length);
		P = new double[y.length];
		X = new double[x.length][];
		for (int i = 0; i < x.length; i++) {
			X[i] = new double[x[i].length];
			System.arraycopy(x[i], 0, X[i], 0, x[i].length);
		}
		B = new double[x[0].length];
	}

	public static void main(String[] args) {
		double[] Y = { 1, 1, 1, 1, 1, 0, 0, 0, 0, 0 };
		double[][] X = { { 1, 15, 4 }, { 1, 30, 14 }, { 1, 31, 16 }, { 1, 31, 11 }, { 1, 32, 17 }, { 1, 29, 10 },
				{ 1, 30, 8 }, { 1, 31, 12 }, { 1, 32, 6 }, { 1, 40, 7 } };
		double[] b = { 0, 0, 0 };
		LogisticRegression LogReg = new LogisticRegression(Y, X);
		LogReg.MLE();
		// RealMatrix Matrix_X = new RealMatrixImpl(X);
		// RealMatrix Matrix_XT = Matrix_X.transpose();
		// int iteration = 20;
		// int iter = 0;
		// RealMatrix B_old = new RealMatrixImpl(b);
		// RealMatrix B_new = new RealMatrixImpl(b);
		// do {
		// B_old = B_new;
		// b = B_old.getColumn(0);
		// RealMatrix Matrix_W = new RealMatrixImpl(getWMatrix(X, b));
		// RealMatrix Matrix_XT_W = Matrix_XT.multiply(Matrix_W);
		// RealMatrix Matrix_XT_W_X = Matrix_XT_W.multiply(Matrix_X);
		// RealMatrix Inv_XT_W_X = Matrix_XT_W_X.inverse();
		// RealMatrix Inv_XT_W_X_XT = Inv_XT_W_X.multiply(Matrix_XT);
		// RealMatrix Vector_Pre = new RealMatrixImpl(predicted(Y, X, b));
		// RealMatrix Vector_H = Inv_XT_W_X_XT.multiply(Vector_Pre);
		// B_new = B_old.add(Vector_H);
		// System.out.println(iter + "--->" + B_new);
		// } while (iter++ < iteration && !stop(B_old, B_new));
		// System.out.println();
	}

	public void MLE() {
		RealMatrix Matrix_X = new RealMatrixImpl(X);
		RealMatrix Matrix_XT = Matrix_X.transpose();
		int iteration = 20;
		int iter = 0;
		RealMatrix B_old = new RealMatrixImpl(B);
		RealMatrix B_new = new RealMatrixImpl(B);
		calculate_P(B);
		LogLikelihood_old = convergence();
		do {
			B_old = B_new;
			LogLikelihood_old = LogLikelihood_new;
			B = B_old.getColumn(0);
			RealMatrix Matrix_W = new RealMatrixImpl(getWMatrix());
			RealMatrix Matrix_XT_W = Matrix_XT.multiply(Matrix_W);
			RealMatrix Matrix_XT_W_X = Matrix_XT_W.multiply(Matrix_X);
			RealMatrix Inv_XT_W_X = Matrix_XT_W_X.inverse();
			RealMatrix Inv_XT_W_X_XT = Inv_XT_W_X.multiply(Matrix_XT);
			RealMatrix Vector_Res = new RealMatrixImpl(Residual());
			RealMatrix Vector_H = Inv_XT_W_X_XT.multiply(Vector_Res);
			B_new = B_old.add(Vector_H);
			calculate_P(B_new.getColumn(0));
			LogLikelihood_new = convergence();
			System.out.println(iter + "--->" + B_new);
			System.out.println(LogLikelihood_old - LogLikelihood_new);
		} while (iter++ < iteration && Math.abs(LogLikelihood_old - LogLikelihood_new) > threshold);
		System.out.println();
	}

	public void calculate_P(double[] b) {
		for (int i = 0; i < X.length; i++) {
			double s = 0;
			for (int j = 0; j < X[i].length ; j++) {
				s += X[i][j] * b[j];
			}
			P[i] = Math.exp(s)/(1 + Math.exp(s));
		}
	}

	public double[] Residual() {
		double[] res = new double[Y.length];
		for (int i = 0; i < P.length; i++) {
			res[i] = Y[i] - P[i];
		}
		return res;
	}

	public double[][] getWMatrix() {
		double[][] W = new double[P.length][P.length];
		for (int i = 0; i < P.length; i++) {
			W[i][i] = P[i] * (1 - P[i]);
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

	public double convergence() {
		boolean converg = true;
		double L = 0;
		for (int i = 0; i < P.length; i++) {
			L += Y[i] * Math.log10(P[i]) + (1 - Y[i]) * Math.log10(1 - P[i]);
		}
		return L * (-1);
	}
}
