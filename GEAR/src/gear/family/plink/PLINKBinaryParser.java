package gear.family.plink;

import gear.family.pedigree.file.BEDReader;
import gear.family.pedigree.file.BIMReader;

public class PLINKBinaryParser extends PLINKParser {

	public static final int HOMOZYGOTE_FIRST = 0x0;
	public static final int HETEROZYGOTE = 0x2;
	public static final int HOMOZYGOTE_SECOND = 0x3;
	public static final int MISSING_GENOTYPE = 0x1;

	protected String famFile;

	public PLINKBinaryParser(String ped, String map, String fam) {
		super(ped, map);
		famFile = fam;
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
	public void parse() {
		mapData = new BIMReader(mapFile);
		ParseMapFile();

		pedData = new BEDReader(famFile, snpFilter.getWorkingSNP().length, mapData);
		pedData.setHeader(false);
		ParsePedFile();
		pedData.cleanup();
	}

	@Override
	public void ParsePedFile() {
		pedData.parseLinkage(pedigreeFile, mapData.getMarkerNumberOriginal(), snpFilter.getWorkingSNP());
	}
}
