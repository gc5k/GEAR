package family.pedigree.file;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import util.NewIt;


public class MapFile {
	private static final String DELIMITER = "\\s+";
	private ArrayList<SNP> snpList;
	private ArrayList<Integer> badline;
	private String mf = null;
	private File mapfile = null;

	public MapFile(String m) {
		mf = m;
	}

	public void parseMap() {
		mapfile = new File(mf);

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
				if (tokens.length < 4) {
					if (badline == null)
						badline = NewIt.newArrayList();
					badline.add(new Integer(c));
				}
				String chr = tokens[0];
				String name = tokens[1];
				// Skip genetic distance field at tokens[2].
				float dis = Float.parseFloat(tokens[2]);
				int pos = Integer.parseInt(tokens[3]);
				snpList.add(new SNP(chr, name, dis, pos));
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
	
	public void setMarker(int l) {
		snpList = NewIt.newArrayList();
		snpList.ensureCapacity(l);
		for(int i = 0; i < l; i++ ) {
			StringBuilder sb = new StringBuilder();
			sb.append("snp");
			sb.append(i);
			snpList.add(new SNP(sb.toString(), "-1", -1f, -1));
		}
	}
	
	public ArrayList<SNP> getMarkerList() {
		return snpList;
	}
	
	public int getMarkerNumber() {
		return snpList.size();
	}
}
