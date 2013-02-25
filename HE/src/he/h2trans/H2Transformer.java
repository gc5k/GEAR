package he.h2trans;
import java.text.DecimalFormat;

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

	public H2Transformer(Parameter p) {
		cal_k = p.cal_k;
		cal_hl = p.cal_hl;
		cal_ho = p.cal_ho;
		cal_cc = p.cal_cc;
		cal_h2_se = p.cal_h2_se;
		P = cal_cc[0]/(cal_cc[0] + cal_cc[1]);
		sb.append("cal_k: " + cal_k + "\n");
		sb.append("cal_hl: " + cal_hl + "\n");
		sb.append("cal_ho: " + cal_ho + "\n");
		sb.append("cal_P: " + P + "\n");
	}

	public void H2() {
		DecimalFormat fmt = new DecimalFormat("#.###E0");

		NormalDistributionImpl Norm = new NormalDistributionImpl();
		double q = 0;
		try {
			q = Norm.inverseCumulativeProbability(1-cal_k);
		} catch (MathException e) {
			e.printStackTrace();
		}
		double z = 1 / (Math.sqrt(2*3.1416926)) * Math.exp(-q*q/2);
		double h2 = 0;
		if(cal_ho!=0) {
			h2 = cal_ho * cal_k * (1-cal_k) * cal_k * (1-cal_k) / (z * z * P * (1-P));
			double hl_se = 0;
			if(cal_h2_se > 0) {
				hl_se = cal_h2_se * (cal_k * (1-cal_k) * cal_k * (1-cal_k))/( z * z * P * (1-P));
				sb.append("h2(l): " + fmt.format(h2) + "\t" + fmt.format(hl_se) + "\n");
			} else {
				sb.append("h2(l): " + fmt.format(h2) + "\n");
			}
		} else {
			h2 = cal_hl * z * z * P * (1-P) / (cal_k * (1- cal_k) * cal_k * (1-cal_k));
			double ho_se = 0;
			if(cal_h2_se > 0) {
				ho_se = cal_h2_se / ( (cal_k * (1-cal_k) * cal_k * (1-cal_k))/( z * z * P * (1-P)) );
				sb.append("h2(o): " + fmt.format(h2) + "\t" + fmt.format(ho_se) + "\n");
			} else {
				sb.append("h2(o): " + fmt.format(h2) + "\n");
			}
		}
		System.out.println(sb);
	}
}
