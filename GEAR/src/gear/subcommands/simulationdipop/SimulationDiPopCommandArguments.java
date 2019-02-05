package gear.subcommands.simulationdipop;

import gear.subcommands.CommandArguments;
import gear.util.FileUtil;
import gear.util.Logger;

public class SimulationDiPopCommandArguments extends CommandArguments {

	public int getRep() {
		return this.rep;
	}

	public void setSampleSize(String[] n) {
		for (int i = 0; i < 2; i++) {
			N[i] = Integer.parseInt(n[i]);
			if (N[i] < 1) {
				Logger.printUserLog("Sample size " + N + " is too small. GEAR quit.");
				System.exit(0);
			}
		}
	}

	public int[] getSampleSize() {
		return N;
	}

	public int getCaseSize() {
		return N[0];
	}

	public int getControlSize() {
		return N[1];
	}

	public void setMarkerNum(String[] m) {
		M = new int[m.length];
		for (int i = 0; i < M.length; i++) {
			M[i] = Integer.parseInt(m[i]);
			if (M[i] < 1) {
				Logger.printUserLog("Marker number " + M[i] + " is too small. GEAR quit.");
				System.exit(0);
			}
		}
	}

	public int[] getMarkerNum() {
		return M;
	}

	public int getTotalMarkerNum() {
		int m = 0; 
		for(int i = 0; i < M.length; i++) {
			m += M[i];
		}
		return m;
	}

	public void setPlainEffect(double e) {
		polyEffect = e;
		isPlainEffect = true;
		isPolyEffect = false;
		isPolyEffectSort = false;
		isPolyEffectFile = false;
	}

	public boolean isPlainEffect() {
		return isPlainEffect;
	}

	public double getPolyEffect() {
		return polyEffect;
	}

	public void setPolyEffect() {
		isPlainEffect = false;
		isPolyEffect = true;
		isPolyEffectSort = false;
		isPolyEffectFile = false;
	}

	public boolean isPolyEffect() {
		return isPolyEffect;
	}

	public void setPolyEffectSort() {
		isPlainEffect = false;
		isPolyEffect = false;
		isPolyEffectSort = true;
		isPolyEffectFile = false;
	}

	public boolean isPolyEffectSort() {
		return isPolyEffectSort;
	}

	public void setPolyEffectFile(String f) {
		FileUtil.exists(f);
		polyEffectFile = f;

		isPlainEffect = false;
		isPolyEffect = false;
		isPolyEffectSort = false;
		isPolyEffectFile = true;
	}

	public boolean isPolyEffectFile() {
		return isPolyEffectFile;
	}

	public String getPolyEffectFile() {
		return polyEffectFile;
	}

	public void setLD(double ld) {
		this.ld = ld;
		if (this.ld < -1 || this.ld > 1) {
			Logger.printUserLog("LD should be between -1 and 1. GEAR quit.");
		}
		isPlainLD = true;
		isRandLD = false;
	}

	public void setLDRange(String[] ld) {
		ldRangeLow = Double.parseDouble(ld[0]);
		ldRangeHigh = Double.parseDouble(ld[1]);
		if (ldRangeLow <= -1 || ldRangeHigh >= 1) {
			Logger.printUserError("LD (Lewontin's) frequency is out of range: " + ldRangeLow + "--" + ldRangeHigh);
		}
		if (ldRangeLow > ldRangeHigh) {
			double t = ldRangeLow;
			ldRangeLow = ldRangeHigh;
			ldRangeLow = t;
		}
		isPlainLD = false;
		isRandLD = true;
	}

	public double getLDRangeLow() {
		return ldRangeLow;
	}

	public double getLDRangeHigh() {
		return ldRangeHigh;
	}

	public boolean isPlainLD() {
		return isPlainLD;
	}

	public double getLD() {
		return ld;
	}

	public void setRandLD() {
		isPlainLD = false;
		isRandLD = true;
	}

	public boolean isRandLD() {
		return isRandLD;
	}

	public void setHsq(double h) {
		hsq = h;
		if (hsq < 0 || hsq > 0.99) {
			Logger.printUserLog("hsq should be between 0 ~ 1. GEAR quit.");
		}
	}

	public double getHsq() {
		return hsq;
	}

	public void setMakeBed() {
		isMakeBed = true;
	}

	public boolean isMakeBed() {
		return isMakeBed;
	}

	public void setFst(String[] fst) {

		if (fst.length != M.length) {
			Logger.printUserLog("The number of Fst does not match the subset of markers. GEAR quit.");
			System.exit(1);
		}

		Fst = new double[fst.length];
		for (int i = 0; i < Fst.length; i++) {
			Fst[i] = Double.parseDouble(fst[i]);
			if (i > 0 && Fst[i] < Fst[i-1]) {
				Logger.printUserLog("Fst should be set ascendingly. GEAR quit.");
				System.exit(1);
			}
		}
		return;
	}

	public double[] getFst() {
		return Fst;
	}

	public void setFreqRange(String[] rf) {
		FreqRangeLow = Double.parseDouble(rf[0]);
		FreqRangeHigh = Double.parseDouble(rf[1]);
		if (FreqRangeLow <= 0 || FreqRangeHigh >= 1) {
			Logger.printUserError("Allele frequency is out of range: " + FreqRangeLow + "--" + FreqRangeHigh);
		}
		if (FreqRangeLow > FreqRangeHigh) {
			double t = FreqRangeLow;
			FreqRangeLow = FreqRangeHigh;
			FreqRangeLow = t;
		}
	}

	public double getFreqRangeLow() {
		return FreqRangeLow;
	}

	public double getFreqRangeHigh() {
		return FreqRangeHigh;
	}

	private int[] N = {100, 100};
	private int[] M = {1000, 100};
	private double[] Fst = {0.01, 0.05};
//	private int nullM = 0;

	private double FreqRangeLow = 0.01;
	private double FreqRangeHigh = 0.5;

	private double polyEffect = 1;
	private boolean isPlainEffect = true;
	private boolean isPolyEffect = false;
	private boolean isPolyEffectSort = false;
	private boolean isPolyEffectFile = false;
	private String polyEffectFile = null;

	private double ld = 0;
	private boolean isPlainLD = true;
	private boolean isRandLD = true;

	private double ldRangeLow = -1;
	private double ldRangeHigh = 1;

	private double hsq = 0.5;
	private boolean isMakeBed = false;

	private int rep = 1;

}
