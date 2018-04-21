package gear.subcommands.ebatchgwas;

import gear.subcommands.CommandArguments;

public class EbatchGWASCommandArguments extends CommandArguments {
	public void setEV(String ev) {
		this.ev = Integer.parseInt(ev);
	}

	public int getEV() {
		return ev;
	}

	public void setDom() {
		isDom = true;
		isEpi = false;
	}

	public boolean isDom() {
		return isDom;
	}

	public void setEpi() {
		isEpi = true;
		isDom = false;
	}

	public boolean isEpi() {
		return isEpi;
	}

	public void setInbred() {
		adjvar = true;
	}

	public boolean isInbred() {
		return adjvar;
	}

	private int ev = 1;
	private boolean isDom = false;
	private boolean isEpi = false;
	private boolean adjvar = false;

}
