package gear.subcommands.propc;

import java.util.ArrayList;

import gear.subcommands.CommandArguments;
import gear.subcommands.profile.ProfileCommand;
import gear.subcommands.profile.ProfileCommandArguments;
import gear.util.BufferedReader;
import gear.util.FileUtil;
import gear.util.NewIt;

public class ProPCCommandArguments extends CommandArguments
{
	public void setBatch(String batch)
	{
		FileUtil.exists(batch);
		BatchFile = batch;
		
		bFile = NewIt.newArrayList();
		BufferedReader reader = BufferedReader.openTextFile(batch, "proPC Batch");

		String[] tokens = null;
		while((tokens = reader.readTokens())!=null)
		{
			String sbed = tokens[0] + ".bed";
			FileUtil.exists(sbed);
			String sbim = tokens[0] + ".bim";
			FileUtil.exists(sbim);
			String sfam = tokens[0] + ".fam";
			FileUtil.exists(sfam);
			bFile.add(tokens[0]);
			System.out.println(tokens[0]);
		}
	}

	public String getBatchFile()
	{
		return BatchFile;
	}

	public ArrayList<String> getbFile()
	{
		return bFile;
	}

	public void setGreedy(boolean greedy)
	{
		this.isGreedy = greedy;
	}

	public boolean isGreedy()
	{
		return isGreedy;
	}

	public ProfileCommandArguments getProfileCommandArguments()
	{
		return profileCommandArguments;
	}

	public void setProfileCommandArguments(ProfileCommandArguments profileCommandArguments)
	{
		this.profileCommandArguments = profileCommandArguments;
	}

	public boolean isGreedy = false;
	public String BatchFile = null;
	public ArrayList<String> bFile = NewIt.newArrayList();
	
	private ProfileCommandArguments profileCommandArguments;
}
