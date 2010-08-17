package LogP;

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
	private double degree_freedom;
	private String key;

	public PointMappingStatistic(int c, int i, int w, double l, 
			double a,double a_sd, double a_t, double a_t_p, 
			double d, double d_sd, double d_t, double d_t_p, double df) {
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
		degree_freedom = df;
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

	public double get_tStatistic_dominant() {
		return D_t;
	}

	public double get_LOD() {
		return lod;
	}
}
