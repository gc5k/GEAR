package gear.family.plink;

import gear.family.pedigree.file.MapFile;
import gear.family.pedigree.file.PedigreeFile;
import gear.family.qc.colqc.SNPFilter;
import gear.subcommands.CommandArguments;
import gear.util.Logger;

public abstract class PLINKParser {
	public static PLINKParser create(CommandArguments cmdArgs) {
		PLINKParser pp = null;
		if (cmdArgs.isbFile()) {
			pp = new PLINKBinaryParser(cmdArgs.getBed(), cmdArgs.getBim(), cmdArgs.getFam());
		} else if (cmdArgs.getFile() != null) {
			pp = new PLINKTextParser(cmdArgs.getPed(), cmdArgs.getMap());
		} else {
			Logger.printUserError("PLINK format file is not specified");
			System.exit(1);
		}
		pp.cmdArgs = cmdArgs;
		return pp;
	}

	public static PLINKParser parse(CommandArguments cmdArgs) {
		PLINKParser pp = create(cmdArgs);
		if (pp != null)
			pp.parse();
		return pp;
	}
	
	protected CommandArguments cmdArgs;

	protected MapFile mapData = null;
	protected PedigreeFile pedigreeData = null;
	protected SNPFilter snpFilter;

	protected PLINKParser(PedigreeFile pedigreeData, MapFile mapData) {
		this.pedigreeData = pedigreeData;
		this.mapData = mapData;
		pedigreeData.setHeader(mapData.getFilename() == null);
		snpFilter = new SNPFilter(mapData);
	}
	
	public void parseSmallFiles() {
		if (mapData.getFilename() != null)
			parseMapFile();
		if (cmdArgs != null) {
			snpFilter.filter(cmdArgs);
			mapData.setWorkingSNPs(snpFilter.getWorkingSNPs());
		}
	}

	public void parse() {
		parseSmallFiles();
		parsePedigreeFile();
	}

	protected void parseMapFile() {
		Logger.printUserLog("Reading '%s'.", mapData.getFilename());
		mapData.parseMap();
		Logger.printUserLog("Read %d marker(s) from '%s'",
				mapData.getMarkerNumberOriginal(), mapData.getFilename());
		Logger.printUserLog("Selected marker(s): " + mapData.getMarkerNumber());
	}

	public void parsePedigreeFile() {
		Logger.printUserLog("Reading '" + pedigreeData.getFilename() + "'.");
		pedigreeData.parseLinkage(snpFilter.getWorkingSNPs());
		Logger.printUserLog("Individual number: " + pedigreeData.getNumIndividuals());
	}

	public PedigreeFile getPedigreeData() {
		return pedigreeData;
	}

	public MapFile getMapData() {
		return mapData;
	}

	public SNPFilter getSNPFilter() {
		return snpFilter;
	}
}
