package gear.subcommands.oath.pick;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;

import gear.subcommands.CommandArguments;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;

public class OathPickCommandArguments extends CommandArguments {

	public void setPickList(String pfile) {
		FileUtil.exists(pfile);
		pickF = pfile;
		BufferedReader br = FileUtil.FileOpen(pickF);
		String oathF = null;
		try {
			while ((oathF = br.readLine()) != null) {
				FileUtil.exists(oathF);
				pickList.add(oathF);
			}
		} catch (IOException e) {
			Logger.handleException(e, "Could not find '" + oathF + "'.");
		}
	}

	public ArrayList<String> getPickList() {
		return pickList;
	}

	public String getPickFile() {
		return pickF;
	}

	private String pickF = null;
	private ArrayList<String> pickList = NewIt.newArrayList();
}
