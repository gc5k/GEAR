package gear.subcommands.grmstat;

import gear.subcommands.CommandArguments;

public class GRMStatCommandArguments extends CommandArguments {

	public String getGrmBin() {
		return grmBin;
	}

	public void setGrmBin(String grmBin) {
		this.grmBin = grmBin;
	}

	public String getGrmText() {
		return grmText;
	}

	public void setGrmText(String grmText) {
		this.grmText = grmText;
	}

	public String getGrmGZ() {
		return grmGZ;
	}

	public void setGrmGZ(String grmGZ) {
		this.grmGZ = grmGZ;
	}

	public String getGrmID() {
		return grmID;
	}

	public void setGrmID(String grmID) {
		this.grmID = grmID;
	}

	private String grmBin;
	private String grmText; // root name of the GRM files
	private String grmGZ; // root name of the GRM files
	private String grmID;

}
