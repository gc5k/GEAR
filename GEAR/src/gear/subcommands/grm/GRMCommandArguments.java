package gear.subcommands.grm;

import java.nio.file.Files;

import gear.subcommands.CommandArguments;
import gear.util.FileUtil;

public class GRMCommandArguments extends CommandArguments {
	private boolean isGZ = true;
	private boolean isVar = false;
	private boolean isDom = false;
	private boolean isInbred = false;
	private boolean isInbredList = false;
	private String InbredFile = null;

	public void setGZ() {
		isGZ = true;
	}

	public void setTxt() {
		isGZ = false;
	}

	public boolean isGZ() {
		return isGZ;
	}

	public void setAdjVar() {
		isVar = true;
	}

	public boolean isAdjVar() {
		return isVar;
	}

	public void setDom() {
		isDom = true;
	}

	public boolean isDom() {
		return isDom;
	}

	public void setInbred() {
		isInbred = true;
	}

	public boolean isInbred() {
		return isInbred;
	}

	public void setInbedList(String inbredL) {
		FileUtil.exists(inbredL);
		isInbredList = true;
		InbredFile = inbredL;
	}
	
	public boolean isInbredList() {
		return isInbredList;
	}
	
	public String getInbredList() {
		return InbredFile;
	}
}
