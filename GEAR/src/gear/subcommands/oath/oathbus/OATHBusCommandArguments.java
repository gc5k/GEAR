package gear.subcommands.oath.oathbus;

import gear.subcommands.CommandArguments;
import gear.util.FileUtil;
import gear.util.Logger;

public class OATHBusCommandArguments extends CommandArguments {
	private String covFile = null;
	private int[] covIdx = { 0 };
	private String removeInter = "11";

	public void setCovFile(String cFile) {
		FileUtil.exists(cFile);
		covFile = cFile;
	}

	public String getCovFile() {
		return covFile;
	}

	public void setCovNumber(String[] cIdx) {
		covIdx = new int[cIdx.length];
		for (int i = 0; i < covIdx.length; i++) {
			covIdx[i] = Integer.parseInt(cIdx[i]);
			if (covIdx[i] < 1) {
				Logger.printUserLog(covIdx[i] + "< 1. Incorrect index for covar-number.");
				Logger.printUserLog("GEAR quittted.");
				System.exit(1);
			}
			covIdx[i]--;
		}
	}

	public int[] getCovNumber() {
		return covIdx;
	}

	public void setRemoveInter(String str) {
		removeInter = str;
		if (removeInter.length() != 2)
			System.exit(1);
	}

	public String RemoveInter() {
		return removeInter;
	}
}
