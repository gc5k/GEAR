package gear.subcommands.fpc;

import java.util.ArrayList;

import gear.subcommands.CommandArguments;
import gear.util.BufferedReader;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;

public class FPCCommandArguments extends CommandArguments
{
	public void setFstFile(String fstF)
	{
		FileUtil.exists(fstF);
		fstFile = fstF;
	}

	public String getFstFile()
	{
		return fstFile;
	}

	public void setReference(String[] fstF)
	{
		ref[0] = Integer.parseInt(fstF[0]) - 1;
		ref[1] = Integer.parseInt(fstF[1]) - 1;
		ref[2] = Integer.parseInt(fstF[2]) - 1;
	}

	public int[] getReference()
	{
		return ref;
	}

	public void setCoordinates(String cf)
	{
		FileUtil.exists(cf);
		BufferedReader reader = BufferedReader.openTextFile(cf, "Summary Statistic file");

		int cnt = 0;
		String[] tokens = null;

		while( (tokens = reader.readTokens(2)) != null && cnt <=2 )
		{
			coordinates[cnt][0] = Double.parseDouble(tokens[0]);
			coordinates[cnt][1] = Double.parseDouble(tokens[1]);
			cnt++;
		}

		if (cnt < 2) 
		{
			Logger.printUserLog("Three pairs of coordinates should be provided.");
		}
	}

	public double[][] getCoordinates()
	{
		return coordinates;
	}

	private String fstFile = null;
	private int[] ref = {0, 1, 2};
	private double[][] coordinates = {{0, -2}, {-1 * Math.sqrt(3), 1}, {Math.sqrt(3), 1}};
}
