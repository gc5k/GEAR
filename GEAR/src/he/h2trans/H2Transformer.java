package he.h2trans;

import gear.util.Logger;
import parameter.Parameter;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistributionImpl;


public class H2Transformer {

	private double cal_k;
	private double cal_hl;
	private double cal_ho;
	private double[] cal_cc = {0,0};
	private double cal_h2_se;
	private double P = 1;
	StringBuffer sb = new StringBuffer();

	public H2Transformer() {
		cal_k = Parameter.INSTANCE.cal_k;
		cal_hl = Parameter.INSTANCE.cal_hl;
		cal_ho = Parameter.INSTANCE.cal_ho;
		cal_cc = Parameter.INSTANCE.cal_cc;
		cal_h2_se = Parameter.INSTANCE.cal_h2_se;
		P = cal_cc[0]/(cal_cc[0] + cal_cc[1]);
		sb.append("K (prevalence): " + cal_k + "\n");
		if(Parameter.INSTANCE.cal_hlFlag) {
			sb.append("heritability on the liability scale: " + cal_hl + "\n");
		}
		if(Parameter.INSTANCE.cal_hoFlag) {
			sb.append("heritability on the observed scale: " + cal_ho + "\n");
		}
		sb.append("Proportion of cases: " + P + "\n");
	}

	public void H2() {
//		DecimalFormat fmt = new DecimalFormat("#.###E0");

		NormalDistributionImpl Norm = new NormalDistributionImpl();
		double q = 0;
		try {
			q = Norm.inverseCumulativeProbability(1-cal_k);
		} catch (MathException e) {
			e.printStackTrace();
		}
		double z = 1 / (Math.sqrt(2*3.1416926)) * Math.exp(-q*q/2);
		double h2 = 0;
		if(Parameter.INSTANCE.cal_hoFlag) {
			h2 = cal_ho * cal_k * (1-cal_k) * cal_k * (1-cal_k) / (z * z * P * (1-P));
			double hl_se = 0;
			if(Parameter.INSTANCE.cal_h2seFlag ) {
				hl_se = cal_h2_se * (cal_k * (1-cal_k) * cal_k * (1-cal_k))/( z * z * P * (1-P));
				sb.append("\nThe transformed h2(l): " + h2 + " se: " + hl_se + "\n");
			} else {
				sb.append("\nh2(l): " + h2 + "\n");
			}
		} else if (Parameter.INSTANCE.cal_hlFlag) {
			h2 = cal_hl * z * z * P * (1-P) / (cal_k * (1- cal_k) * cal_k * (1-cal_k));
			double ho_se = 0;
			if(Parameter.INSTANCE.cal_h2seFlag ) {
				ho_se = cal_h2_se / ( (cal_k * (1-cal_k) * cal_k * (1-cal_k))/( z * z * P * (1-P)) );
				sb.append("\nh2(o): " + h2 + " se: " + ho_se + "\n");
			} else {
				sb.append("\nh2(o): " + h2 + "\n");
			}
		}
		Logger.printUserLog(sb.toString());
	}
}
