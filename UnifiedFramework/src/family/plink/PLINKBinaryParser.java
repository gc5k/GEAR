package family.plink;

import java.io.File;
import java.io.IOException;

import family.pedigree.file.BEDReader;
import family.pedigree.file.BIMReader;
import family.pedigree.file.GMDRPhenoFile;

public class PLINKBinaryParser extends PLINKParser {

	protected String FamFile;
	public PLINKBinaryParser(String ped, String phe, String map, String fam) {
		super(ped, phe, map);
		FamFile = fam;

	}

	@Override
	public void initial() {
		mapData = new BIMReader(mapFile);
		phenoData = new GMDRPhenoFile();
		if (mapFile != null) {
			ParseMapFile();
			pedData = new BEDReader(FamFile, mapData.getMarkerNumber(), mapData);
			pedData.setHeader(false);
			ParsePedFile();
		}
		mapData.setPolymorphism(pedData.getAlleleFrequency());
		pedData.cleanup();
		if (phenotypeFile != null) {
			ParsePhenoFile();
		}
	}

	@Override
	public void ParsePedFile() {
		File PedFile = new File(pedigreeFile);
		try {
			pedData.parseLinkage(PedFile, mapData.getMarkerNumber());
		} catch (IOException e) {
			System.err.println("Pedgree file initialization exception.");
			e.printStackTrace(System.err);
		}
	}
}
