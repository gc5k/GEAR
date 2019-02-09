package gear.family.plink;

import gear.family.pedigree.file.BEDReader;
import gear.family.pedigree.file.BIMReader;

public class PLINKBinaryParser extends PLINKParser {
	public static final byte HOMOZYGOTE_FIRST = 0b00;
	public static final byte HETEROZYGOTE = 0b10;
	public static final byte HOMOZYGOTE_SECOND = 0b11;
	public static final byte MISSING_GENOTYPE = 0b01;

	protected String famFile;
	private BEDReader bedReader;
	
	private PLINKBinaryParser(String bed, BIMReader bimReader, String fam) {
		super(new BEDReader(bed, fam, bimReader), bimReader);
		bedReader = (BEDReader)pedigreeData;
	}

	public PLINKBinaryParser(String bed, String map, String fam) {
		this(bed, new BIMReader(map), fam);
		famFile = fam;
	}
	
	@Override
	public void parseSmallFiles() {
		super.parseSmallFiles();
		bedReader.prepareToParseGenotypes();
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
}
