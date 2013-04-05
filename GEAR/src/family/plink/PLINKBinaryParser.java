package family.plink;

import java.io.IOException;
import java.util.logging.Level;

import family.pedigree.file.BEDReader;
import family.pedigree.file.BIMReader;
import gear.util.Logger;

public class PLINKBinaryParser extends PLINKParser {
	
	public static final int HOMOZYGOTE_FIRST = 0x0;
	public static final int HETEROZYGOTE = 0x2;
	public static final int HOMOZYGOTE_SECOND = 0x3;
	public static final int MISSING_GENOTYPE = 0x1;

	protected String FamFile;
	
	public PLINKBinaryParser(String ped, String map, String fam) {
		super(ped, map);
		FamFile = fam;
	}
	
	public static int convertToGearGenotype(int plinkGenotype) {
		switch (plinkGenotype) {
		case PLINKBinaryParser.HOMOZYGOTE_FIRST:
			return gear.ConstValues.BINARY_HOMOZYGOTE_FIRST;
		case PLINKBinaryParser.HETEROZYGOTE:
			return gear.ConstValues.BINARY_HETEROZYGOTE;
		case PLINKBinaryParser.HOMOZYGOTE_SECOND:
			return gear.ConstValues.BINARY_HOMOZYGOTYE_SECOND;
		case PLINKBinaryParser.MISSING_GENOTYPE:
			return gear.ConstValues.BINARY_MISSING_GENOTYPE;
		default:
			assert false;
			return gear.ConstValues.BINARY_MISSING_GENOTYPE;
		}
	}

	@Override
	public void Parse() {
		mapData = new BIMReader(mapFile);
		if (mapFile != null) {
			ParseMapFile();
			Logger.printUserLog("Reading '" + mapFile + "'.");
			Logger.printUserLog("Marker Number: " + mapData.getMarkerNumberOriginal());
			pedData = new BEDReader(FamFile, snpFilter.getWorkingSNP().length, mapData);
			pedData.setHeader(false);
			ParsePedFile();
			Logger.printUserLog("Reading '" + pedigreeFile + "'.");
			Logger.printUserLog("Individual Number: " + pedData.getNumIndividuals());
		}
		pedData.cleanup();
	}

	@Override
	public void ParsePedFile() {
		try {
			pedData.parseLinkage(pedigreeFile, mapData.getMarkerNumberOriginal(), snpFilter.getWorkingSNP());
		} catch (IOException e) {
			Logger.printUserError("An exception occurred when parsing the pedigree files.");
			Logger.printUserError("Exception Message: " + e.getMessage());
			Logger.getDevLogger().log(Level.SEVERE, "Parsing pedigree files", e);
			System.exit(1);
		}
	}
}
