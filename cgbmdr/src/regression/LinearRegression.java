package regression;

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
    double sd_residual;
    double SSTO;
    double SSR;
    double SSE;
    double F_statistic;
    double P_F_statistic;
    int df;
    int df_numerator;
    int df_denominator;
    double[] Var;
    double[] residuals;

    public LinearRegression(double[][] x, double [][] y) {
        Predictor = new RealMatrixImpl(x);
        Response = new RealMatrixImpl(y);
        df = Response.getRowDimension() - Predictor.getColumnDimension();
        Var = new double[Predictor.getColumnDimension()];
        residuals = new double[y.length];
        df_numerator = Predictor.getColumnDimension() - 1;
        df_denominator = df;
    }

    public void MLE() {
        RealMatrix X_t = Predictor.transpose();
        RealMatrix X_tX = X_t.multiply(Predictor);
        RealMatrix X_tX_Ivt = X_tX.inverse();
        RealMatrix X_tX_Ivt_Xt = X_tX_Ivt.multiply(X_t);
        estimate = X_tX_Ivt_Xt.multiply(Response);

        RealMatrix fit = Predictor.multiply(estimate);
        RealMatrix residual = Response.subtract(fit);
        residuals = residual.getColumn(0);
        RealMatrix residual_t = residual.transpose();
        RealMatrix sse = residual_t.multiply(residual);
        double d = (Response.getRowDimension()-Predictor.getColumnDimension());
        double mse = sse.getEntry(0, 0)/d;
        RealMatrix FI = X_tX_Ivt;
        TDistribution t = new TDistributionImpl(d);
        for(int i = 0; i < Var.length; i++) {
            Var[i] = mse*FI.getEntry(i, i);
            double t_val = Math.abs(estimate.getEntry(i, 0))/Math.sqrt(Var[i]);
            try {
                double p = t.cumulativeProbability(t_val);
            } catch (MathException E) {
                E.printStackTrace(System.err);
            }
        }

        SummaryStatistics SS_ssto = new SummaryStatisticsImpl();
        SummaryStatistics SS_sse =  new SummaryStatisticsImpl();
        for(int i = 0; i < Response.getRowDimension(); i++) {
            SS_ssto.addValue(Response.getEntry(i, 0));
            SS_sse.addValue(residual.getEntry(i, 0));
        }
        SSTO = SS_ssto.getVariance()*SS_ssto.getN();
        SSE = SS_sse.getVariance()*SS_sse.getN();
        SSR = SSTO - SSE;
        if(df_numerator <= 0) {
            F_statistic = 0;
            P_F_statistic = 1;
        } else {
            F_statistic = (SSR/df_numerator)/(SSE/df_denominator);
            FDistribution fd = new FDistributionImpl(df_numerator, df_denominator);
            try {
                P_F_statistic = 1 - fd.cumulativeProbability(F_statistic);
            } catch (MathException E) {
                E.printStackTrace(System.err);
            }
        }
    }

    public double get_F_Statistic () {
        return F_statistic;
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

    public RealMatrix getEstimate() {
        return estimate;
    }

    public double[] getResiduals() {
        return residuals;
    }
    
    public static void main (String[] args) {
        double[][] data = { {1, 1, 2, 3},
                            {1, 1.3, 2.1, 3.4},
                            {1, 1.2, 1.0, 3.9},
                            {1, 0.8, 2.3, 3.4},
                            {1, 0.9, 2.5, 3.1},
                            {1, 0.98, 3.2, 3.9},
                          };
        double[][] y_data = {{5}, {7}, {6}, {5.6}, {6.7}, {6.5}};
        LinearRegression lm = new LinearRegression(data, y_data);
        lm.MLE();
        double[] res = lm.getResiduals();
        for(int i = 0; i < res.length; i++) {
            System.out.println(res[i]);
        }
    }
}
