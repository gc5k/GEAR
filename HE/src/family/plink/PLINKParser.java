package family.plink;

import java.io.File;
import java.io.IOException;

import test.Test;

import parameter.Parameter;

import family.qc.colqc.SNPFilter;
import family.qc.colqc.SNPFilterI;
import family.qc.colqc.SNPFilterInterface;
import family.pedigree.file.MapFile;
import family.pedigree.file.PedigreeFile;

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
		pedData.setHeader(Parameter.header);

		if (mapFile != null) {//bim
			ParseMapFile();
			Test.LOG.append("reading " + mapFile + ".\n");
			Test.LOG.append(mapData.getMarkerNumberOriginal() + " markers in " + mapFile + ".\n");
			Test.LOG.append(mapData.getMarkerNumber() + " selected markers.\n");
			System.err.println("reading " + mapFile + ".");
			System.err.println(mapData.getMarkerNumberOriginal() + " markers in " + mapFile + ".");
			System.err.println(mapData.getMarkerNumber() + " selected markers.");
			pedData.setHeader(false);
			ParsePedFile();
			Test.LOG.append("reading " + pedigreeFile + ".");
			Test.LOG.append(pedData.getNumIndividuals() + " individuals.\n");
			System.err.println("reading " + pedigreeFile + ".");
			System.err.println(pedData.getNumIndividuals() + " individuals.");
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
		if (Parameter.transFlag) {
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
	 * @param Ped
	 *            the name of the pedigree file
	 * @throws IOException
	 */
	public void ParsePedFile() {

		try {
			pedData.parseLinkage(pedigreeFile, mapData.getMarkerNumberOriginal(), snpFilter.getWorkingSNP());
		} catch (IOException e) {
			System.err.println("Pedgree file initialization exception.");
			e.printStackTrace(System.err);
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
