package gear.qc.snpqc;

import java.text.DecimalFormat;

import gear.subcommands.CommandArguments;
import gear.util.Logger;

public class SNPFilterPostQC {
	private boolean isMAF = false;
	private int mafCnt = 0;
	private boolean isGeno = false;
	private int genoCnt = 0;
	private double genoM = 0.0;
	private int gCnt = 0;

	private boolean isMaxMAF = false;
	private int maxmafCnt = 0;
	private boolean isMAFRange = false;
	private int mafRangeCnt = 0;

	private CommandArguments qcArgs = null;

	public SNPFilterPostQC(CommandArguments cmdArgs) {
		qcArgs = cmdArgs;
		if (qcArgs.isGENO()) {
			isGeno = true;
		}
		if (qcArgs.isMAF()) {
			isMAF = true;
		}
		if (qcArgs.isMaxMAF()) {
			isMaxMAF = true;
		}
		if (qcArgs.isMAFRange()) {
			isMAFRange = true;
		}
	}
	
	public boolean isPassPostQC(double freq, double mRate) {
		double maf = freq < 0.5 ? freq:(1-freq);

		if (isGeno && mRate > qcArgs.getGENO()) {
			genoCnt++;
			return false;
		}
		if (isMAF && maf < qcArgs.getMAF()) {
			mafCnt++;
			return false;
		}
		if (isMaxMAF && maf > qcArgs.getMaxMAF()) {
			maxmafCnt++;
			return false;
		}
		if (isMAFRange) {
			double[][] mafR = qcArgs.getMAFRange();
			for (int j = 0; j < mafR.length; j++) {
				if (maf >= mafR[j][0] && maf <= mafR[j][1]) {
					mafRangeCnt++;
					return false;
				}
			}
		}
		genoM += mRate;
		gCnt++;
		return true;
	}
	
	public void printPostQCSummary() {
		int Cnt = mafCnt + maxmafCnt + mafRangeCnt + genoCnt;
		DecimalFormat fmt1 = new DecimalFormat("0.0000");

		if (isMAF || isMaxMAF || isMAFRange || isGeno) {
			Logger.printUserLog("In total " + Cnt + " SNPs were excluded in QC(s).");			
		}
		if (isGeno) {
			Logger.printUserLog(genoCnt + " SNPs were excluded because did not pass missing rate threshold " + qcArgs.getGENO() + ".");
			Cnt += genoCnt;
		}
		if (isMAF) {
			Logger.printUserLog(mafCnt + " SNPs were excluded because did not pass maf threshold " + qcArgs.getMAF() + ".");
			Cnt += mafCnt;
		}
		if (isMaxMAF) {
			Logger.printUserLog(maxmafCnt + " SNPs were excluded because did not pass max-maf threshod " + qcArgs.getMaxMAF() + ".");
			Cnt += maxmafCnt;
		}
		if (isMAFRange) {
			Logger.printUserLog(mafRangeCnt + " SNPs were excluded because did not pass maf-range thresholds.");
			Cnt += mafRangeCnt;
		}
		Logger.printUserLog("Average genotyping rate is " + fmt1.format(1-genoM/gCnt) + " for " +gCnt + " loci passed QC.");
	}
}
