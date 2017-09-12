package gear.subcommands.ebatchgwas;

import gear.subcommands.CommandArguments;

public class EbatchGWASArguments extends CommandArguments
{
	public void setEV(String ev)
	{
		this.ev = Integer.parseInt(ev);
	}
	
	public int getEV()
	{
		return ev;
	}

	public void setDom(boolean isdom)
	{
		isDom = isdom;
	}

	public boolean isDom()
	{
		return isDom;
	}

	private int ev = 1;
	private boolean isDom = false;

}
