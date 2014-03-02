package gear.subcommands.metawatchdog.encrypt;

import gear.subcommands.CommandArguments;
import gear.subcommands.profile.ProfileCommandArguments;

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

	public ProfileCommandArguments getProfileCommandArguments()
	{
		return profileCommandArguments;
	}

	public void setProfileCommandArguments(ProfileCommandArguments profileCommandArguments)
	{
		this.profileCommandArguments = profileCommandArguments;
	}

	private String mapFile;
	private String encodeFile;
	private ProfileCommandArguments profileCommandArguments;
}
