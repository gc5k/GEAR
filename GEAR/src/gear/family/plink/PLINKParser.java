package gear.family.plink;

import gear.family.pedigree.file.MapFile;
import gear.family.pedigree.file.PedigreeFile;
import gear.family.qc.colqc.SNPFilter;
import gear.subcommands.CommandArguments;
import gear.util.Logger;

public class PLINKParser {
	public static PLINKParser parse(CommandArguments cmdArgs) {
		PLINKParser pp = null;
		if (cmdArgs.isbFile()) {
			pp = new PLINKBinaryParser(cmdArgs.getBed(), cmdArgs.getBim(), cmdArgs.getFam());
		} else if (cmdArgs.getFile() != null) {
			pp = new PLINKParser(cmdArgs.getPed(), cmdArgs.getMap());
		} else {
			return null;
		}
		pp.doParse(cmdArgs);
		return pp;
	}

	protected MapFile mapData = null;
	protected PedigreeFile pedData = null;
	protected SNPFilter snpFilter;
	protected String pedigreeFile;
	protected String mapFile;

	public PLINKParser(String ped, String map) {
		pedigreeFile = ped;
		mapFile = map;
	}

	public void Parse() {
		mapData = new MapFile(mapFile);

		pedData = new PedigreeFile();
		pedData.setHeader(false);

		if (mapFile != null) {// bim
			ParseMapFile();
			Logger.printUserLog("Reading '" + mapFile + "'.");
			Logger.printUserLog("Marker number: " + mapData.getMarkerNumberOriginal());
			Logger.printUserLog("Selected marker number: " + mapData.getMarkerNumber());
			pedData.setHeader(false);
			ParsePedFile();
			Logger.printUserLog("Reading '" + pedigreeFile + "'.");
			Logger.printUserLog("Individual number: " + pedData.getNumIndividuals());
		} else {
			pedData.setHeader(true);
			ParsePedFile();
			mapData.setMarker(pedData.getNumMarker());
		}
		mapData.setPolymorphismMarker(pedData.getPolymorphism());
		pedData.cleanup();

	}

	protected void doParse(CommandArguments cmdArgs) {
		mapData = new MapFile(mapFile);

		pedData = new PedigreeFile();
		pedData.setHeader(false);

		if (mapFile != null) {// bim
			ParseMapFile(cmdArgs);
			Logger.printUserLog("Reading '" + mapFile + "'.");
			Logger.printUserLog("Marker number: " + mapData.getMarkerNumberOriginal());
			Logger.printUserLog("Selected marker number: " + mapData.getMarkerNumber());
			pedData.setHeader(false);
			ParsePedFile();
			Logger.printUserLog("Reading '" + pedigreeFile + "'.");
			Logger.printUserLog("Individual number: " + pedData.getNumIndividuals());
		} else {
			pedData.setHeader(true);
			ParsePedFile();
			mapData.setMarker(pedData.getNumMarker());
		}
		mapData.setPolymorphismMarker(pedData.getPolymorphism());
		pedData.cleanup();

	}

	public void ParseMapFile() {
		mapData.parseMap();
		Logger.printUserLog("Read " + mapData.getMarkerNumberOriginal() + " SNPs from '" + mapFile + "'.");
		snpFilter = new SNPFilter(mapData);
		snpFilter.SelectSNP();
		int[] WSNP = snpFilter.getWorkingSNP();
		mapData.setWSNP(WSNP);
	}

	public void ParseMapFile(CommandArguments cmdArgs) {
		mapData.parseMap();
		Logger.printUserLog("Read " + mapData.getMarkerNumberOriginal() + " SNPs from '" + cmdArgs.getBim() + "'.");
		snpFilter = new SNPFilter(mapData);
		snpFilter.SelectSNP(cmdArgs);
		int[] WSNP = snpFilter.getWorkingSNP();
		mapData.setWSNP(WSNP);
	}

	/**
	 * Initialize basic implementation of the genotype file.
	 * 
	 * @param Ped
	 *            the name of the pedigree file
	 */
	public void ParsePedFile() {
		pedData.parseLinkage(pedigreeFile, mapData.getMarkerNumberOriginal(), snpFilter.getWorkingSNP());
	}

	public PedigreeFile getPedigreeData() {
		return pedData;
	}

	public MapFile getMapData() {
		return mapData;
	}

	public SNPFilter getSNPFilter() {
		return snpFilter;
	}

	public void setAlleleFrequency(double[][] freq) {
		mapData.setAlleleFrequency(freq);
	}
}
