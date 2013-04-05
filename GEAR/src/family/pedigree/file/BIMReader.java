package family.pedigree.file;

import gear.util.Logger;
import gear.util.NewIt;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.logging.Level;

public class BIMReader extends MapFile {

	public BIMReader(String m) {
		super(m);
	}
	
	public void parseMap() {
		mapfile = new File(mf);
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(mapfile));
		} catch (IOException e) {
			Logger.printUserError("Can't open the map file.");
			Logger.printUserError("Exception Message: " + e.getMessage());
			Logger.getDevLogger().log(Level.SEVERE, "Creating BufferedReader", e);
			System.exit(1);
		}

		String line = null;
		try {
			int c = 0;
			while ((line = reader.readLine()) != null) {
				c++;
				String[] tokens = line.split(DELIMITER);
				if (tokens.length < 6) {
					if (badline == null)
						badline = NewIt.newArrayList();
					badline.add(new Integer(c));
				}
				String chr = tokens[0];
				String name = tokens[1];
				// Skip genetic distance field at tokens[2].
				float dis = Float.parseFloat(tokens[2]);
				int pos = Integer.parseInt(tokens[3]);
				addSNP(chr, name, dis, pos, tokens[4].charAt(0), tokens[5].charAt(0));
			}
			reader.close();
		} catch (IOException e) {
			Logger.printUserError("An exception occurred when reading the map file '" + mf + "'.");
			Logger.printUserError("Exception Message: " + e.getMessage());
			Logger.getDevLogger().log(Level.SEVERE, "Parsing the map file", e);
			System.exit(1);
		}

		if (badline != null) {
			Logger.printUserError("problems with the lines below:");
			String badlines = "";
			for (Integer i: badline) {
				badlines += i + ",";
			}
			Logger.printUserError(badlines);
			System.exit(1);
		}
		numMarkerOriginal = snpList.size();
	}

}
