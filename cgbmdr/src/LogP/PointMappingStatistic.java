package LogP;

import org.apache.commons.math.distribution.FDistribution;
import org.apache.commons.math.distribution.FDistributionImpl;

public class PointMappingStatistic {
	private int chr;
	private int interval;
	private int walk;
	private double lod;
	private double A;
	private double A_sd;
	private double A_t;
	private double A_t_p_value_cu;
	private double D;
	private double D_sd;
	private double D_t;
	private double D_t_p_value_cu;
	private double degree_freedom_t;
	private double wald;
	private double p_wald_cu;
	private double degree_freedom_wald;
	private double mse;
	private String key;

	public PointMappingStatistic(int c, int i, int w, double l, 
			double a,double a_sd, double a_t, double a_t_p, 
			double d, double d_sd, double d_t, double d_t_p, double df_t, double wa, double pw, double df_w, double m) {
		chr = c;
		interval = i;
		walk = w;
		lod = l;
		A = a;
		A_sd = a_sd;
		A_t = a_t;
		A_t_p_value_cu = a_t_p;
		D = d;
		D_sd = d_sd;
		D_t = d_t;
		D_t_p_value_cu = d_t_p;
		degree_freedom_t = df_t;
		wald = wa;
		p_wald_cu = pw;
		degree_freedom_wald = df_w;
		mse = m;
		int[] substring = { chr, interval, walk };
		key = Utils.makeKey(substring);
	}

	public String getKey() {
		return key;
	}

	public double get_logP_additive() {
		double logp = 0;
		if (A_t < 0) {
			logp = -1 * Math.log10(A_t_p_value_cu*2);
		} else {
			logp = -1 * Math.log10((1 - A_t_p_value_cu)*2);
		}
		return logp;
	}

	public double get_P_additive() {
		if (A_t < 0) {
			return A_t_p_value_cu * 2;
		} else {
			return (1-A_t_p_value_cu) * 2;
		}
	}

	public double get_logP_dominance() {
		double logp = 0;
		if (D_t < 0) {
			logp = -1 * Math.log10(D_t_p_value_cu*2);
		} else {
			logp = -1 * Math.log10((1 - D_t_p_value_cu)*2);
		}
		return logp;
	}

	public double get_P_dominance() {
		if (A_t < 0) {
			return A_t_p_value_cu * 2;
		} else {
			return (1-A_t_p_value_cu) * 2;
		}
	}

	
	public double get_tStatistic_additive() {
		return A_t;
	}

	public double get_tStatistic_dominance() {
		return D_t;
	}

	public double get_LOD() {
		return lod;
	}
	
	public double get_wald() {
		return wald;
	}

	public double get_P_wald() {
		return 1 - p_wald_cu;
	}
	
	public double get_logP_wald() {
		return -1 * Math.log10(1 - p_wald_cu);
	}
	
	public double get_P_F() {
		double pf = 0;
		double s = (wald/degree_freedom_wald)/mse;
		FDistribution fd = new FDistributionImpl(degree_freedom_wald, degree_freedom_t);
		try {
			pf = 1 - fd.cumulativeProbability(s);
		} catch (Exception E) {
			E.printStackTrace(System.err);
		}
		return pf;
	}
	
	public double get_logP_F() {
		double pf = 0;
		double s = (wald/degree_freedom_wald)/mse;
		FDistribution fd = new FDistributionImpl(degree_freedom_wald, degree_freedom_t);
		try {
			pf = 1 - fd.cumulativeProbability(s);
		} catch (Exception E) {
			E.printStackTrace(System.err);
		}
		return (-1) * Math.log10(pf);
	}
}
