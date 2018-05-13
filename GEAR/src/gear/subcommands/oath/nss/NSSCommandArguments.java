package gear.subcommands.oath.nss;

import gear.subcommands.CommandArguments;
import gear.subcommands.CommandSplitter;
import gear.util.FileUtil;

public class NSSCommandArguments extends CommandArguments {
	private String covFile = null;
	private int[] covIdx = { 0 };

	public void setCovFile(String cFile) {
		FileUtil.exists(cFile);
		covFile = cFile;
	}

	public String getCovFile() {
		return covFile;
	}

	public void setCovarIndex(String[] cIdx) {
		CommandSplitter CS = new CommandSplitter.Builder(cIdx).OPT(NSSCommand.OPT_COVAR_NUMBER).IntMin(1).create();
		covIdx = CS.ParseToInt();
	}

	public void setCovNumber(int[] cIdx) {
		covIdx = new int[cIdx.length];
		System.arraycopy(cIdx, 0, covIdx, 0, cIdx.length);
	}

	public int[] getCovNumber() {
		return covIdx;
	}
}
