package gear.subcommands.ibd;

import gear.subcommands.CommandArguments;

public class IBDCommandArguments extends CommandArguments
{
	private boolean isGZ = true;
	private boolean scale = false;

	public void setGZ()
	{
		isGZ = true;
	}

	public void setTxt()
	{
		isGZ = false;
	}
	
	public boolean isGZ()
	{
		return isGZ;
	}

	public void setScale() 
	{
		scale = true;
	}
	
	public boolean isScale()
	{
		return scale;
	}

}
