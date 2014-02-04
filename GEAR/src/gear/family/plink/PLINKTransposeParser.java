package gear.family.plink;

import gear.family.pedigree.file.MapFile;
import gear.family.pedigree.file.TransposePedigreeReader;

public class PLINKTransposeParser extends PLINKParser
{

	protected String FamFile;

	public PLINKTransposeParser(String ped, String Fam, String map)
	{
		super(ped, null);
		FamFile = Fam;
	}

	public void Parse()
	{
		mapData = new MapFile(null);
		pedData = new TransposePedigreeReader(pedigreeFile, FamFile, mapData);
		pedData.setHeader(false);

		pedData.setHeader(false);
		ParsePedFile();

		mapData.setPolymorphismMarker(pedData.getPolymorphism());
		pedData.cleanup();
	}

	public void ParsePedFile()
	{
		pedData.parseLinkage(pedigreeFile, 0, snpFilter.getWorkingSNP());
	}
}
