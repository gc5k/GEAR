package family.plink;

import java.io.IOException;
import test.Test;
import family.pedigree.file.BEDReader;
import family.pedigree.file.BIMReader;

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

	@Override
	public void Parse() {
		mapData = new BIMReader(mapFile);
		if (mapFile != null) {
			ParseMapFile();
			Test.LOG.append("reading " + mapFile + ".\n");
			Test.LOG.append(mapData.getMarkerNumberOriginal() + " markers in " + mapFile + ".\n");
			System.err.println("reading " + mapFile + ".");
			System.err.println(mapData.getMarkerNumberOriginal() + " markers in " + mapFile + ".");
			pedData = new BEDReader(FamFile, snpFilter.getWorkingSNP().length, mapData);
			pedData.setHeader(false);
			ParsePedFile();
			Test.LOG.append("reading " + pedigreeFile + ".");
			Test.LOG.append(pedData.getNumIndividuals() + " individuals.\n");
			System.err.println("reading " + pedigreeFile + ".");
			System.err.println(pedData.getNumIndividuals() + " individuals.");
		}
		pedData.cleanup();
	}

	@Override
	public void ParsePedFile() {
		try {
			pedData.parseLinkage(pedigreeFile, mapData.getMarkerNumberOriginal(), snpFilter.getWorkingSNP());
		} catch (IOException e) {
			System.err.println("Pedgree file initialization exception.");
			e.printStackTrace(System.err);
		}
	}
}
