package family.plink;

import java.io.File;
import java.io.IOException;

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
		}
		if (mapFile != null) {
			ParseMapFile();
			pedData.setHeader(false);
			ParsePedFile();
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
