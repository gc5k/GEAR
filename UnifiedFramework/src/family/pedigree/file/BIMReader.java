package family.pedigree.file;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import util.NewIt;

public class BIMReader extends MapFile {

	public BIMReader(String m) {
		super(m);
	}
	
	public void parseMap() {
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(mapfile));
		} catch (IOException E) {
			System.err.println("can't open map file\n");
		}
		snpList = NewIt.newArrayList();

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
				snpList.add(new SNP(chr, name, dis, pos, tokens[4].charAt(0), tokens[5].charAt(0)));
			}
			reader.close();
		} catch (IOException E) {
			System.err.println("bad map file");
		}

		if(badline != null) {
			System.err.println( "problems with the lines below:");
			for(Integer i: badline) {
				System.err.println( i + ",");
			}
			System.exit(0);
		}
	}

}
