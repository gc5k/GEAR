package family.plink;

import java.io.IOException;
import java.util.logging.Level;

import parameter.Parameter;

import family.qc.colqc.SNPFilter;
import family.qc.colqc.SNPFilterI;
import family.qc.colqc.SNPFilterInterface;
import family.pedigree.file.MapFile;
import family.pedigree.file.PedigreeFile;
import gear.util.Logger;

public class PLINKParser {

	protected MapFile mapData = null;
	protected PedigreeFile pedData = null;
//	protected PhenotypeFile phenoData = null;
	protected SNPFilterInterface snpFilter;
	protected String pedigreeFile;
	protected String phenotypeFile;
	protected String mapFile;

	public PLINKParser(String ped, String map) {
		pedigreeFile = ped;
		mapFile = map;
	}

	public void Parse() {
		mapData = new MapFile(mapFile);

		pedData = new PedigreeFile();
		pedData.setHeader(false);

		if (mapFile != null) {//bim
			ParseMapFile();
			Logger.printUserLog("Reading '" + mapFile + "'.");
			Logger.printUserLog("Marker Numeber: " + mapData.getMarkerNumberOriginal());
			Logger.printUserLog("Selected Marker Number: " + mapData.getMarkerNumber());
			pedData.setHeader(false);
			ParsePedFile();
			Logger.printUserLog("Reading '" + pedigreeFile + "'.");
			Logger.printUserLog("Individual Number: " + pedData.getNumIndividuals());
		} else {
			pedData.setHeader(true);
			ParsePedFile();
			mapData.setMarker(pedData.getNumMarker());
		}
		mapData.setPolymorphismMarker(pedData.getPolymorphism());
		pedData.cleanup();

	}

	public void ParseMapFile() {
		if (mapFile != null) {
			mapData.parseMap();
		}
		if (Parameter.INSTANCE.transFlag) {
			snpFilter = new SNPFilterI(mapData);
		} else {
			snpFilter = new SNPFilter(mapData);
		}
		snpFilter.Select();
		int[] WSNP = snpFilter.getWorkingSNP();
		mapData.setWSNP(WSNP);
	}

	/**
	 * Initialize basic implementation of the genotype file.
	 * 
	 * @param Ped	the name of the pedigree file
	 */
	public void ParsePedFile() {
		try {
			pedData.parseLinkage(pedigreeFile, mapData.getMarkerNumberOriginal(), snpFilter.getWorkingSNP());
		} catch (IOException e) {
			Logger.printUserError("An exception occurred when parsing the pedigree files.");
			Logger.printUserError("Exception message: " + e.getMessage());
			Logger.getDevLogger().log(Level.SEVERE, "Parsing pedigree files", e);
			System.exit(1);
		}
	}

	public PedigreeFile getPedigreeData() {
		return pedData;
	}

	public MapFile getMapData() {
		return mapData;
	}

	public SNPFilterInterface getSNPFilter() {
		return snpFilter;
	}

	public void setAlleleFrequency(double[][] freq) {
		mapData.setAlleleFrequency(freq);
	}
}
