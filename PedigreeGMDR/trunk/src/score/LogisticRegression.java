package score;
import java.util.*;
import org.apache.commons.math.linear.*;

public class LogisticRegression {
	double[] Y;
	double[][] X;
	double[] B;
	double threshold = 0.001;
	double LogLikelihood_old;
	double LogLikelihood_new;
	LogisticRegression(double[]y, double[][] x) {
		Y = new double[y.length];
		System.arraycopy(y, 0, Y, 0, Y.length);
		X = new double[x.length][];
		for (int i = 0; i < x.length; i++) {
			X[i] = new double[x[i].length];
			System.arraycopy(x[i], 0, X[i], 0, X[i].length);
		}
		B = new double[X[0].length];
	}

	public static void main(String[] args) {
		double[] Y = {  1, 1, 1, 1, 1, 0, 0, 0, 0, 0};
		double[][] X = { { 1, 15, 4}, { 1, 30, 14}, { 1, 31, 16}, { 1, 31, 11}, { 1, 32, 17}, 
				         { 1, 29, 10}, { 1, 30, 8}, { 1, 31, 12}, { 1, 32, 6}, { 1, 40, 7} };
		double[] b = { 0, 0, 0};
		LogisticRegression LogReg = new LogisticRegression(Y, X);
		LogReg.MLE();
//		RealMatrix Matrix_X = new RealMatrixImpl(X);
//		RealMatrix Matrix_XT = Matrix_X.transpose();
//		int iteration = 20;
//		int iter = 0;
//		RealMatrix B_old = new RealMatrixImpl(b);
//		RealMatrix B_new = new RealMatrixImpl(b);
//		do {
//			B_old = B_new;
//			b = B_old.getColumn(0);
//			RealMatrix Matrix_W = new RealMatrixImpl(getWMatrix(X, b));
//			RealMatrix Matrix_XT_W = Matrix_XT.multiply(Matrix_W);
//			RealMatrix Matrix_XT_W_X = Matrix_XT_W.multiply(Matrix_X);
//			RealMatrix Inv_XT_W_X = Matrix_XT_W_X.inverse();
//			RealMatrix Inv_XT_W_X_XT = Inv_XT_W_X.multiply(Matrix_XT);
//			RealMatrix Vector_Pre = new RealMatrixImpl(predicted(Y, X, b));
//			RealMatrix Vector_H = Inv_XT_W_X_XT.multiply(Vector_Pre);
//			B_new = B_old.add(Vector_H);
//			System.out.println(iter + "--->" + B_new);
//		} while (iter++ < iteration && !stop(B_old, B_new));
//		System.out.println();
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
			RealMatrix Matrix_W = new RealMatrixImpl(getWMatrix(X, B));
			RealMatrix Matrix_XT_W = Matrix_XT.multiply(Matrix_W);
			RealMatrix Matrix_XT_W_X = Matrix_XT_W.multiply(Matrix_X);
			RealMatrix Inv_XT_W_X = Matrix_XT_W_X.inverse();
			RealMatrix Inv_XT_W_X_XT = Inv_XT_W_X.multiply(Matrix_XT);
			RealMatrix Vector_Pre = new RealMatrixImpl(predicted(Y, X, B));
			RealMatrix Vector_H = Inv_XT_W_X_XT.multiply(Vector_Pre);
			B_new = B_old.add(Vector_H);
			LogLikelihood_new = convergence(B_new);
			System.out.println(iter + "--->" + B_new);
			System.out.println(LogLikelihood_old - LogLikelihood_new);
		} while(iter++ < iteration && Math.abs(LogLikelihood_old - LogLikelihood_new) > threshold);
		System.out.println();
	}

	public double[] predicted(double[] y, double[][] x, double[] b) {
		double[] res = new double[y.length];
		for (int i = 0; i < x.length; i++) {
			double s = 0;
			for (int j = 0; j < x[i].length; j++) {
				s += x[i][j] * b[j];
			}
			res[i] = y[i] - Math.exp(s) / (Math.exp(s) + 1);
		}
		return res;
	}

	public double[] _predicted(double[][] x, double[] b) {
		double[] y = new double[x.length];
		for (int i = 0; i < x.length; i++) {
			double s = 0;
			for (int j = 0; j < x[i].length; j++) {
				s += x[i][j] * b[j];
			}
			y[i] = 1 / (Math.exp(s) + 1);
		}
		return y;
	}

	public double[][] getWMatrix(double[][] x, double[] b) {
		double[][] W = new double[x.length][x.length];
		for (int i = 0; i < x.length; i++) {
			double s = 0;
			for (int j = 0; j < x[i].length; j++) {
				s += x[i][j] * b[j];
			}
			double p = Math.exp(s) / (1 + Math.exp(s));
			W[i][i] = p * (1-p);
		}
		return W;
	}
	
	public boolean stop(RealMatrix b_old, RealMatrix b_new) {
		boolean stopit = true;
		for(int i = 0; i < b_old.getRowDimension(); i++) {
			if(Math.abs(b_old.getEntry(i, 0) - b_new.getEntry(i, 0))>threshold) {
				stopit = false;
				break;
			}
		}
		return stopit;
	}
	
	public double convergence(RealMatrix b) {
		boolean converg = true;
		double L = 0;
		for(int i = 0; i < X.length; i++) {
			double s = 0;
			for(int j = 0; j < X[i].length; j++) {
				s += X[i][j] * b.getEntry(j, 0);
			}
			double p = Math.exp(s) / (1 + Math.exp(s));
			L += Y[i] * Math.log10(p) + (1-Y[i]) * Math.log10(1-p);
		}
		return L * (-1);
	}
}
