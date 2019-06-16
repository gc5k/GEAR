package gear.subcommands.wgrmA;

import gear.subcommands.CommandArguments;
import gear.util.FileUtil;

public class WGRMACommandArguments extends CommandArguments {
	private boolean isGZ = true;
	private boolean isVar = false;
	private boolean isDom = false;
	private boolean isDomOnly = false;
	private boolean isInbred = false;
	private boolean isInbredList = false;
	private String InbredFile = null;
	private boolean isMemLow = false;
	private boolean isGUI = false;

	private boolean isVanRaden = false;
	private boolean isWeight = false;
	private String wFile = null;

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

	public void setGUI() {
		isGUI = true;
	}

	public boolean isGUI() {
		return isGUI;
	}

	public void setVanRaden() {
		isVanRaden = true;
		isWeight = false;
	}

	public boolean isWeight() {
		return isWeight;
	}

	public boolean isVanRaden() {
		return isVanRaden;
	}

	public void setWeightFile(String wF) {
		FileUtil.exists(wF);
		wFile = wF;
		isWeight = true;
	}

	public String getWeightFile() {
		return wFile;
	}

	public void setDomOnly() {
		isDomOnly = true;
	}
	
	public boolean isDomOnly() {
		return isDomOnly;
	}

	public void setMemLow() {
		isMemLow = true;	
	}
	
	public boolean isMemLow() {
		return isMemLow;
	}
}
