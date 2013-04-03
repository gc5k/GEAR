package gear.util.stat;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistributionImpl;

public class Z {

	public static double OddsRatioTest(double f1, double f2, double n1, double n2) {
		double odds = Math.log((f1*(1-f2))/((1-f1)*f2));
		double se = Math.sqrt(1/(f1*n1) + 1/((1-f1)*n1) + 1/(f2*n2) + 1/((1-f2)*n2));
		double T = odds/se;
		return T;
	}

	public static double OddsRatioTestPvalueOneTail(double f1, double f2, double n1, double n2) {
		double T = OddsRatioTest(f1, f2, n1, n2);

		NormalDistributionImpl nd = new NormalDistributionImpl();
		double p = 1;
		try {
			p = 1 - nd.cumulativeProbability(T);
		} catch (MathException e) {
			e.printStackTrace();
		}

		return p;
	}

	public static double OddsRatioTestPvalueTwoTail(double f1, double f2, double n1, double n2) {
		double T = OddsRatioTest(f1, f2, n1, n2);
		NormalDistributionImpl nd = new NormalDistributionImpl();
		double p = 1;
		try {
			p = 1 - nd.cumulativeProbability(T);
		} catch (MathException e) {
			e.printStackTrace();
		}
		
		if (2 * p > 1) {
			p = 2*(1-p);
		} else {
			p = 2*p;
		}

		return p;
	}

}
