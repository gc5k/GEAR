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

	private int ev = 1;

}
