package he;

import parameter.Parameter;

public class Lambda {
	protected double L = 0;
	protected double XYProd = 0;
	protected double SX = 0;
	protected double SY = 0;
	protected double EX = 0;
	protected double EY = 0;
	protected double XXProd = 0;
	protected double N;

	protected double cov;
	protected double var;

	private int mode = -2; //-2 for squared difference; -1 for cross product
	
	public Lambda () {
		if (Parameter.INSTANCE.heType[Parameter.INSTANCE.he_sd]) {
			mode = -2;
		} else if (Parameter.INSTANCE.heType[Parameter.INSTANCE.he_ss]) {
			mode = -2;
		} else if (Parameter.INSTANCE.heType[Parameter.INSTANCE.he_cp]) {
			mode = -1;
		}
	}

	public void calLambda() {
		
		EX = SX/N;
		EY = SY/N;
		cov = XYProd/N - EX*EY;
		var = XXProd/N - EX*EX;

		L = cov/var/mode;
		
	}
	
	public double getCov() {
		return cov;
	}

	public double getVar() {
		return var;
	}
	
	public double getLambda(double twice_vg) {
		L = cov/var/twice_vg;
		return L;
	}
	
}
