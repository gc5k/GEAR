package gear.subcommands;

import java.util.HashSet;

import gear.util.Logger;
import gear.util.NewIt;

public abstract class CommandArguments {
	public void setFile(String file) {
		this.file = file;
	}

	public String getFile() {
		return file;
	}

	public String getPed() {
		return getFile() + ".ped";
	}

	public String getMap() {
		return getFile() + ".map";
	}

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

	public String getOutRoot() {
		return outRoot;
	}

	public void setOutRoot(String outRoot) {
		this.outRoot = outRoot;
	}

	public long getSeed() {
		return seed;
	}

	public void setSeed(long seed) {
		this.seed = seed;
	}

	public void setKeepFile(String kF) {
		this.keepFile = kF;
		isKeepFile = true;
	}

	public String getKeepFile() {
		return keepFile;
	}

	public boolean isKeepFile() {
		return isKeepFile;
	}

	public void setRemoveFile(String rF) {
		this.removeFile = rF;
		isRemoveFile = true;
	}

	public String getRemoveFile() {
		return removeFile;
	}
	
	public boolean isRemoveFile() {
		return isRemoveFile;
	}

	public void setExtractFile(String eF) {
		this.extractFile = eF;
		isExtractFile = true;
	}

	public String getExtractFile() {
		return extractFile;
	}

	public boolean isExtractFile() {
		return isExtractFile;
	}

	public void setExcludeFile(String exF) {
		this.excludeFile = exF;
		isExcludeFile = true;
	}

	public String getExcludeFile() {
		return excludeFile;
	}

	public boolean isExcludeFile() {
		return isExcludeFile;
	}

	public void setChr(String[] c) {
		chr = NewIt.newHashSet();
		for (int i = 0; i < c.length; i++) {
			chr.add(c[i]);
		}
		isChr = true;
	}

	public HashSet<String> getChr() {
		return chr;
	}

	public boolean isChr() {
		return isChr;
	}

	public void setNotChr(String[] c) {
		not_chr = NewIt.newHashSet();
		for (int i = 0; i < c.length; i++) {
			not_chr.add(c[i]);
		}
		isNotChr = true;
	}

	public HashSet<String> getNotChr() {
		return not_chr;
	}

	public boolean isNotChr() {
		return isNotChr;
	}

	public void setMAF(String m) {
		maf = Double.parseDouble(m);
		if (maf < 0 || maf > 0.5) {
			Logger.printUserLog("incorrect --maf " + m + ". It should be > 0 and < 0.5");
		}
		isMAF = true;
	}

	public double getMAF() {
		return maf;
	}

	public boolean isMAF() {
		return isMAF;
	}

	public void setMaxMAF(String m) {
		max_maf = Double.parseDouble(m);
		if (max_maf < 0 || max_maf > 0.5) {
			Logger.printUserLog("incorrect --max-maf " + m + ". It should be > 0 and < 0.5");
		}
		isMaxMAF = true;
	}

	public double getMaxMAF() {
		return max_maf;
	}

	public boolean isMaxMAF() {
		return isMaxMAF;
	}

	public void setGENO(String m) {
		geno = Double.parseDouble(m);
		if (geno < 0 || geno > 1) {
			Logger.printUserLog("incorrect --geno " + m + ". It should be > 0 and < 1");
		}
		isGENO = true;
	}

	public double getGENO() {
		return geno;
	}

	public boolean isGENO() {
		return isGENO;
	}

	private String file;
	private String bfile;
	private String outRoot;

	private String keepFile = null;
	private boolean isKeepFile = false;

	private String removeFile = null;
	private boolean isRemoveFile = false;

	private String extractFile = null;
	private boolean isExtractFile = false;

	private String excludeFile = null;
	private boolean isExcludeFile = false;

	private HashSet<String> chr = null;
	private boolean isChr = false;
	private HashSet<String> not_chr = null;
	private boolean isNotChr = false;

	private double maf = 0.01;
	private boolean isMAF = false;
	private double max_maf = 0.51;
	private boolean isMaxMAF = false;
	private double geno = 0.1;
	private boolean isGENO = false;
	private long seed;
}
