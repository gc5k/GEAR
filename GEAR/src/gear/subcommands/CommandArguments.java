package gear.subcommands;

import java.util.HashSet;

import gear.util.Logger;
import gear.util.NewIt;

public abstract class CommandArguments {
	public void setFile(String file) {
		this.file = file;
		this.isFile = true;
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

	public boolean isFile() {
		return isFile;
	}

	public void setBFile(String bfile) {
		this.bfile = bfile;
		this.isbFile = true;
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

	public boolean isbFile() {
		return isbFile;
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

	public void setKeepFamFile(String kFamF) {
		this.keepFamFile = kFamF;
		isKeepFamFile = true;
	}

	public String getKeepFamFile() {
		return keepFamFile;
	}

	public boolean isKeepFamFile() {
		return isKeepFamFile;
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

	public void setRemoveFamFile(String rFamF) {
		this.removeFamFile = rFamF;
		isRemoveFamFile = true;
	}

	public String getRemoveFamFile() {
		return removeFamFile;
	}
	
	public boolean isRemoveFamFile() {
		return isRemoveFamFile;
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
		CommandSplitter CS = new CommandSplitter.Builder(c).OPT("--chr").create();
		chr.addAll(CS.ParseToStr());
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
		CommandSplitter CS = new CommandSplitter.Builder(c).OPT("--not-chr").create();
		not_chr.addAll(CS.ParseToStr());
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
			Logger.printUserLog("incorrect --maf " + m + ". It should be > 0 and < 0.5.");
			System.exit(0);
		}
		isMAF = true;
	}

	public double getMAF() {
		return maf;
	}

	public boolean isMAF() {
		return isMAF;
	}

	public void setMAFRange(String[] mafRg) {
		CommandSplitter CS = new CommandSplitter.Builder(mafRg).DoubleMin(0).DoubleMax(0.5).OPT(Command.OPT_MAF_RANGE_LONG).create();
		maf_range = CS.ParseToDouble();
		isMafRange = true;
	}

	public void setMAFRange2(double[][] mafRg2) {
		maf_range = mafRg2;
		isMafRange = true;
	}

	public boolean isMAFRange() {
		return isMafRange;
	}
	
	public double[][] getMAFRange() {
		return maf_range;
	}

	public void setMaxMAF(String m) {
		max_maf = Double.parseDouble(m);
		if (max_maf < 0 || max_maf > 0.5) {
			Logger.printUserLog("incorrect --max-maf " + m + ". It should be > 0 and < 0.5.");
			System.exit(0);
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
			Logger.printUserLog("incorrect --geno " + m + ". It should be > 0 and < 1.");
			System.exit(0);
		}
		isGENO = true;
	}

	public double getGENO() {
		return geno;
	}

	public boolean isGENO() {
		return isGENO;
	}

	public void setZeroVar() {
		isZeroVar = true;
	}

	public boolean isZeroVar() {
		return isZeroVar;
	}

	public String getPhenotypeFile() {
		return this.pheFile;
	}

	public void setPhenotypeFile(String pheFile) {
		this.pheFile = pheFile;
	}

	public void setPhenotypeIndex(String[] pheIndex) {
		CommandSplitter CS = new CommandSplitter.Builder(pheIndex).OPT(Command.OPT_MPHE).IntMin(1).create();
		this.mPheno = CS.ParseToInt();
		for(int i = 0; i < this.mPheno.length; i++) {
			this.mPheno[i]--;
			if (this.mPheno[i] < 0) {
				Logger.printUserLog("incorrect --mphe " + (mPheno[i]+1) + ". It should be > 0.");
				System.exit(0);
			}
		}
	}

	public void setPhenotypeIndex(int pheIndex) {
		this.mPheno[0] = pheIndex - 1;
		if (this.mPheno[0] < 0) {
			Logger.printUserLog("incorrect --mphe " + (mPheno[0]+1) + ". It should be > 0.");
			System.exit(0);
		}
	}

	public int getSelectedPhenotype(int i) {
		return this.mPheno[i];
	}

	public int[] getSelectedPhenotype() {
		return this.mPheno;
	}

	public void setThreadNum(String threadNum) {
		thread_num = Integer.parseInt(threadNum);
		if (thread_num < 1) {
			Logger.printUserLog("incorrect --thread-num " + thread_num + ". It should be > 0.");
			System.exit(0);
		} else {
			int cpuTotal = Runtime.getRuntime().availableProcessors();
			if (thread_num > cpuTotal) {
				Logger.printUserLog("Only " + cpuTotal + " cpus are available. Thread number is set to " + cpuTotal + ".");
				thread_num = cpuTotal;
			}
		}
		isThread_num = true;
	}

	public void setThreadGreedy(String threadGreedy) {

		int cpuTotal = Runtime.getRuntime().availableProcessors();
		int thread_G= Integer.parseInt(threadGreedy);
		if (thread_G > 0) {
			Logger.printUserLog("incorrect --thread-greedy " + thread_G + ". It should be < 0.");
			System.exit(0);
		} else {
			if ((cpuTotal-thread_G) <1) {
				Logger.printUserLog("At least 1 cpu should be used. Thread number is set to " + 1 + ".");
				thread_num = 1;
			} else {
				thread_num = cpuTotal - thread_G;
			}
		}
		isThread_num = true;
	}

	public boolean isThreadNum() {
		return isThread_num;
	}

	public int getThreadNum() {
		return thread_num;
	}

	private String outRoot;

	private String pheFile;
	private int[] mPheno = { 0 };	
	
	private String file;
	private boolean isFile = false;

	private String bfile;
	private boolean isbFile = false;


	private String keepFamFile = null;
	private boolean isKeepFamFile = false;

	private String removeFamFile = null;
	private boolean isRemoveFamFile = false;

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

	private double[][] maf_range;
	private boolean isMafRange = false;

	private double geno = 0.1;
	private boolean isGENO = false;

	private boolean isZeroVar = false;
	private long seed;
	
	private boolean isThread_num = false;
	private int thread_num = 1;
}
