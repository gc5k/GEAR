package gear.subcommands.simuqtreal;

import gear.subcommands.CommandArguments;
import gear.util.FileUtil;
import gear.util.Logger;

public class SimulationQTRealCommandArguments extends CommandArguments {

	public String getBFile() {
		return bfile;
	}

	public void setBFile(String bfile) {
		this.bfile = bfile;
	}

	public String getBed() {
		return getBFile() + ".bed";
	}

	public String getBim() {
		return getBFile() + ".bim";
	}

	public String getFam() {
		return getBFile() + ".fam";
	}

	public void setRep(String rep) {
		this.rep = Integer.parseInt(rep);
		if (this.rep < 1) {
			Logger.printUserLog("Replication should be greater than 0. GEAR quit.");
			System.exit(0);
		}
	}

	public int getRep() {
		return this.rep;
	}

	public void setPlainEffect(double e) {
		polyEffect = e;
		isPlainEffect = true;
		isPolyEffect = false;
		isPolyEffectSort = false;
		isPolyEffectFile = false;
		isRefEffectFile = false;
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
		isRefEffectFile = false;
	}

	public boolean isPolyEffect() {
		return isPolyEffect;
	}

	public void setPolyEffectSort() {
		isPlainEffect = false;
		isPolyEffect = false;
		isPolyEffectSort = true;
		isPolyEffectFile = false;
		isRefEffectFile = false;
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
		isRefEffectFile = false;
	}

	public boolean isPolyEffectFile() {
		return isPolyEffectFile;
	}

	public String getPolyEffectFile() {
		return polyEffectFile;
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

	public void setRefEffectFile(String pFile) {
		FileUtil.exists(pFile);
		polyRefEffectFile = pFile;
		isPlainEffect = false;
		isPolyEffect = false;
		isPolyEffectSort = false;
		isPolyEffectFile = false;
		isRefEffectFile = true;
	}

	public String getRefEffectFile() {
		return polyRefEffectFile;
	}

	public boolean isRefEffectFile() {
		return isRefEffectFile;
	}
	// private int N = 100;
	// private int M = 100;
	// private int nullM = 0;

	private double polyEffect = 1;
	private boolean isPlainEffect = false;
	private boolean isPolyEffect = true;
	private boolean isPolyEffectSort = false;
	private boolean isPolyEffectFile = false;
	private boolean isRefEffectFile = false;
	private String polyEffectFile = null;
	private String polyRefEffectFile = null;

	private double hsq = 0.5;
	private boolean isMakeBed = false;

	private int rep = 1;
	private String bfile;

}
