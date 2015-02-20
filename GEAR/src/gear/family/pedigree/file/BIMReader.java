package gear.family.pedigree.file;

import gear.util.Logger;
import gear.util.NewIt;
import gear.util.BufferedReader;

//import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashSet;

public class BIMReader extends MapFile
{

	public BIMReader(String m)
	{
		super(m);
	}

	public void parseMap()
	{
		mapfile = new File(mf);
		BufferedReader reader = BufferedReader.openTextFile(mf, "bim");
		String[] tokens = null;
		int c = 0;
		HashSet<String> markerSet = NewIt.newHashSet();
		while ((tokens = reader.readTokens(6)) != null)
		{
			c++;
			if (tokens.length < 6)
			{
				String badline = new String();
				for(int i = 0; i < tokens.length; i++)
				{
					badline += tokens[i] + "  ";
				}
				Logger.printUserLog("Line " + c + " has fewer than 6 elements. '" + badline + "'");
				Logger.printUserError("Line " + c + " has fewer than 6 elements. '" + badline + "'");
				Logger.printUserLog("GEAR quitted.");
				System.exit(0);
			}
			String chr = tokens[0];
			String name = tokens[1];
			float dis = Float.parseFloat(tokens[2]);
			int pos = Integer.parseInt(tokens[3]);
			addSNP(chr, name, dis, pos, tokens[4].charAt(0),
					tokens[5].charAt(0));
			if (markerSet.contains(name))
			{
				Logger.printUserLog("'"+ name + "' duplicated.");
			}
			else
			{
				markerSet.add(name);
			}
		} while ((tokens = reader.readTokens(6)) != null);

		numMarkerOriginal = snpList.size();
	}

}
