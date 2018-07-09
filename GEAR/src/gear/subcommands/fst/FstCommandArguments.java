package gear.subcommands.fst;

import gear.subcommands.CommandArguments;
import gear.util.FileUtil;

public class FSTCommandArguments extends CommandArguments {

	
	public void setGroup(String gFile) {
		FileUtil.exists(gFile);
		groupFile = gFile;
	}

	public String getGroupFile() {
		return groupFile;
	}

	private String groupFile = null;
}
