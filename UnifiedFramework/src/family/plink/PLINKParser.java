package family.plink;

import java.io.File;
import java.io.IOException;

import test.Test;

import admixture.parameter.Parameter;

import family.mdr.filter.SNPFilter;
import family.mdr.filter.SNPFilterI;
import family.mdr.filter.SNPFilterInterface;
import family.pedigree.file.MapFile;
import family.pedigree.file.PedigreeFile;
import family.pedigree.file.PhenotypeFile;

public class PLINKParser {

	protected MapFile mapData = null;
	protected PedigreeFile pedData = null;
	protected PhenotypeFile phenoData = null;
	protected SNPFilterInterface snpFilter;
	protected String pedigreeFile;
	protected String phenotypeFile;
	protected String mapFile;

	public PLINKParser(String ped, String map, String phe) {
		pedigreeFile = ped;
		mapFile = map;
		phenotypeFile = phe;
	}

	public void Parse() {
		mapData = new MapFile(mapFile);

		pedData = new PedigreeFile();
		pedData.setHeader(Parameter.header);

		if (phenotypeFile != null) {
			phenoData = new PhenotypeFile();
			ParsePhenoFile();
			Test.LOG.append("reading " + phenotypeFile + ".\n");
			System.err.println("reading " + phenotypeFile + ".");
			Test.LOG.append(phenoData.getNumTraits() + " traits in " + phenotypeFile + ".\n");
			System.err.println(phenoData.getNumTraits() + " traits in " + phenotypeFile + ".");
		}
		if (mapFile != null) {
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
//		mapData.setPolymorphism(pedData.getPolymorphism(), pedData.getAlleleFrequency());
		mapData.setPolymorphismMarker(pedData.getPolymorphism());
		pedData.cleanup();

	}

	public void ParseMapFile() {
		if (mapFile != null) {
			mapData.parseMap();
		}
		if (Parameter.x) {
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

	/**
	 * Initialize basic implementation of the phenotype file.
	 * 
	 * @param Pheno
	 *            the name of the pedigree file
	 * @throws IOException
	 */
	public void ParsePhenoFile() {
		File PheFile = new File(phenotypeFile);
		try {
			phenoData.parsePhenotype(PheFile);
		} catch (IOException e) {
			System.err.println("Pheno file initialization exception.");
			e.printStackTrace(System.err);
		}
	}

	public PedigreeFile getPedigreeData() {
		return pedData;
	}

	public PhenotypeFile getPhenotypeData() {
		return phenoData;
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
