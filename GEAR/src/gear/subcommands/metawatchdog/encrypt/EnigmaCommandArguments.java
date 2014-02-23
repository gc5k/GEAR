package gear.subcommands.metawatchdog.encrypt;

import gear.subcommands.CommandArguments;

public class EnigmaCommandArguments extends CommandArguments
{
	public String getMapFile()
	{
		return mapFile;
	}
	
	public void setMapFile(String mapFile)
	{
		this.mapFile = mapFile;
	}
	
	public String getEncodeFile()
	{
		return encodeFile;
	}

	public void setEncodeFile(String ecodeFile)
	{
		this.encodeFile = ecodeFile;
	}

	private String mapFile;
	private String encodeFile;
}
