package regression;

import im.GenomeScan;
import im.population.IMPopulation;
import java.util.ArrayList;
import org.apache.commons.math.distribution.FDistribution;
import org.apache.commons.math.distribution.FDistributionImpl;
import org.apache.commons.math.distribution.TDistribution;
import org.apache.commons.math.distribution.TDistributionImpl;
import org.apache.commons.math.stat.descriptive.SummaryStatistics;
import org.apache.commons.math.stat.descriptive.SummaryStatisticsImpl;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealMatrixImpl;
import org.apache.commons.math.MathException;

/**
 * 
 * @author Guo-Bo Chen chenguobo@gmail.com
 */
public class LinearRegression {

	RealMatrix Response;
	RealMatrix Predictor;
	RealMatrix estimate;
	RealMatrix residual;
	RealMatrix FisherInformation;
	double SSTO;
	double SSR;
	double SSE;
	double mse;
	double F_statistic;
	double P_F_statistic;
	double logScore;
	int df_Predictor;
	int df_residual;
	double[] p_value_t_test;
	double[] Var;
	double[] tStatic;
	double[] residuals;

	public LinearRegression(double[][] x, double[][] y) {
		Predictor = new RealMatrixImpl(x);
		Response = new RealMatrixImpl(y);
		df_residual = Response.getRowDimension()
				- Predictor.getColumnDimension();
		Var = new double[Predictor.getColumnDimension()];
		p_value_t_test = new double[Predictor.getColumnDimension()];
		tStatic = new double[Predictor.getColumnDimension()];
		residuals = new double[y.length];
		df_Predictor = Predictor.getColumnDimension() - 1;
	}

	public void MLE() {
		RealMatrix X_t = Predictor.transpose();
		RealMatrix X_tX = X_t.multiply(Predictor);
		RealMatrix X_tX_Ivt = X_tX.inverse();
		RealMatrix X_tX_Ivt_Xt = X_tX_Ivt.multiply(X_t);
		estimate = X_tX_Ivt_Xt.multiply(Response);
		RealMatrix fit = Predictor.multiply(estimate);
		residual = Response.subtract(fit);
		residuals = residual.getColumn(0);
		RealMatrix residual_t = residual.transpose();
		RealMatrix sse = residual_t.multiply(residual);
		mse = sse.getEntry(0, 0) / df_residual;

		RealMatrix FI = X_tX_Ivt.scalarMultiply(mse);
		FisherInformation = X_tX.scalarMultiply(1/mse);
		TDistribution t = new TDistributionImpl(df_residual);
		for (int i = 0; i < Var.length; i++) {
			Var[i] = FI.getEntry(i, i);
			tStatic[i] = estimate.getEntry(i, 0) / Math.sqrt(Var[i]);
			try {
				p_value_t_test[i] = t.cumulativeProbability(tStatic[i]);
			} catch (MathException E) {
				E.printStackTrace(System.err);
			}
		}

		SummaryStatistics SS_ssto = new SummaryStatisticsImpl();
		SummaryStatistics SS_sse = new SummaryStatisticsImpl();
		for (int i = 0; i < Response.getRowDimension(); i++) {
			SS_ssto.addValue(Response.getEntry(i, 0));
			SS_sse.addValue(residual.getEntry(i, 0));
		}
		SSTO = SS_ssto.getVariance() * SS_ssto.getN();
		SSE = SS_sse.getVariance() * SS_sse.getN();
		SSR = SSTO - SSE;
		if (df_Predictor <= 0) {
			F_statistic = 0;
			P_F_statistic = 1;
		} else {
			F_statistic = (SSR / df_Predictor) / (SSE / df_residual);
			FDistribution fd = new FDistributionImpl(df_Predictor, df_residual);
			try {
				P_F_statistic = 1 - fd.cumulativeProbability(F_statistic);
			} catch (MathException E) {
				E.printStackTrace(System.err);
			}
		}
	}

	public double[][] quasiResidual(ArrayList selectedMarker, int interval) {
		double[][] Y_res = new double[Response.getRowDimension()][1];
		if (selectedMarker == null) {
			for (int i = 0; i < Y_res.length; i++) {
				Y_res[i][0] = Response.getEntry(i, 0);
				for (int j = 0; j < estimate.getRowDimension(); j++) {
					if (j == (interval + 1) || j == (interval + 2)) {
						continue;
					}
					Y_res[i][0] -= Predictor.getEntry(i, j)
							* estimate.getEntry(j, 0);
				}
			}
		} else {
			for (int i = 0; i < Y_res.length; i++) {
				Y_res[i][0] = Response.getEntry(i, 0);
				for (int j = 0; j < estimate.getRowDimension(); j++) {
					if (j != 0) {
						int mi = ((Integer) selectedMarker.get(j - 1))
								.intValue();
						if (mi == interval || mi == (interval + 1)) {
							continue;
						}
					}
					Y_res[i][0] -= Predictor.getEntry(i, j)
							* estimate.getEntry(j, 0);
				}
			}
		}
		return Y_res;
	}

	public double[][] getResponse() {
		double[][] y = new double[Response.getRowDimension()][Response
				.getColumnDimension()];
		for (int i = 0; i < Response.getRowDimension(); i++) {
			for (int j = 0; j < Response.getColumnDimension(); j++) {
				y[i][j] = Response.getEntry(i, j);
			}
		}
		return y;
	}

	public double getSD(int idx) {
		return Var[idx];
	}

	public double DFResidual() {
		return df_residual;
	}

	public double DFPredictor() {
		return df_Predictor;
	}

	public double DFTotal() {
		return df_residual + df_Predictor + 1;
	}

	public double get_F_Statistic() {
		return F_statistic;
	}

	public double getTStatic(int idx) {
		return tStatic[idx];
	}

	public double getP_F() {
		return P_F_statistic;
	}

	public double getSSTO() {
		return SSTO;
	}

	public double getSSR() {
		return SSR;
	}

	public double getSSE() {
		return SSE;
	}

	public double getMSE() {
		return mse;
	}

	public RealMatrix getEstimate() {
		return estimate;
	}

	public double getCoefficient(int idx) {
		return estimate.getEntry(idx, 0);
	}

	public RealMatrix Y() {
		return Response;
	}

	public RealMatrix X() {
		return Predictor;
	}

	public double[] getResiduals() {
		return residuals;
	}

	public double getWald(int n) {
		double d = 0;
		double[] e = new double[n];
		for (int i = 0; i < e.length; i++) {
			e[i] = estimate.getEntry(1 + i, 0);
		}
		double[] t = new double[n];
		for (int i = 0; i < e.length; i++) {
			for (int j = 0; j < e.length; j++) {
				t[i] += e[j] * FisherInformation.getEntry(1 + j, i + 1);
			}
		}
		for (int i = 0; i < t.length; i++) {
			d += t[i] * e[i];
		}
		return d;
	}

	public double getPValueTTest(int idx) {
		return p_value_t_test[idx];
	}

	public RealMatrix getResidual() {
		return residual;
	}

	public static void main(String[] args) {
		double[][] data = { { 1, 1, 2, 3 }, { 1, 1.3, 2.1, 3.4 },
				{ 1, 1.2, 1.0, 3.9 }, { 1, 0.8, 2.3, 3.4 },
				{ 1, 0.9, 2.5, 3.1 }, { 1, 0.98, 3.2, 3.9 }, };
		double[][] y_data = { { 5 }, { 7 }, { 6 }, { 5.6 }, { 6.7 }, { 6.5 } };
		LinearRegression lm = new LinearRegression(data, y_data);
		lm.MLE();
		double[] res = lm.getResiduals();
		for (int i = 0; i < res.length; i++) {
			System.out.println(res[i]);
		}
	}
}
