package family.plink;

import java.io.File;

import family.pedigree.file.BEDReader;
import family.pedigree.file.BIMReader;

public class PLINKBinaryParser extends PLINKParser {

	protected String FamFile;
	public PLINKBinaryParser(String ped, String phe, String map, String fam) {
		super(ped, phe, map);
		FamFile = fam;
	}

	protected void initial() {
		mapData = new BIMReader(mapFile);

		if (mapFile != null) {
			ParseMapFile();
			pedData = new BEDReader(FamFile, mapData.getMarkerNumber());
			pedData.setHeader(false);
			ParsePedFile();
		}
		mapData.setPolymorphism(pedData.getPolymorphism(), pedData.getAlleleFrequency());
		pedData.cleanup();
		if (phenotypeFile != null) {
			ParsePhenoFile();
		}
	}
	
	
}
