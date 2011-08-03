package family.plink;

import java.io.File;
import java.io.IOException;

import admixture.parameter.Parameter;

import family.pedigree.file.GMDRPhenoFile;
import family.pedigree.file.GMDRPhenoFileException;
import family.pedigree.file.MDRPedFileException;
import family.pedigree.file.MapFile;
import family.pedigree.file.PedigreeFile;

public class PLINKParser {
	protected MapFile mapData;
	protected PedigreeFile pedData;
	protected GMDRPhenoFile phenoData;

	protected String pedigreeFile;
	protected String phenotypeFile;
	protected String mapFile;
	
	public PLINKParser (String ped, String phe, String map) {
		pedigreeFile = ped;
		phenotypeFile = phe;
		mapFile = map;

		mapData = new MapFile(mapFile);
		phenoData = new GMDRPhenoFile();
		pedData = new PedigreeFile();
		pedData.setHeader(Parameter.header);
		initial();
	}

	private void initial() {
		if(mapFile != null) {
			ParseMapFile();
			pedData.setHeader(false);
			ParsePedFile();
		} else {
			pedData.setHeader(true);
			ParsePedFile();
			mapData.setMarker(pedData.getNumMarker());
		}

		if (phenotypeFile != null) {
			ParsePhenoFile();
		}
	}

	public void ParseMapFile() {
		if(mapFile != null) {
			mapData.parseMap();
		}
	}
	
	/**
	 * Initialize basic implementation of the genotype file.
	 * 
	 * @param Ped
	 *            the name of the pedigree file
	 * @throws IOException
	 */
	public void ParsePedFile() {
		File PedFile = new File(pedigreeFile);
		try {
			pedData.Initial(PedFile);
			pedData.parseLinkage();
		} catch (MDRPedFileException e) {
			System.err.println("PedFile initial exception.");
			e.printStackTrace(System.err);
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
			phenoData.Initial(PheFile);
			phenoData.parsePhenotype();
		} catch (GMDRPhenoFileException e) {
			System.err.println("Pheno file initialization exception.");
			e.printStackTrace(System.err);
		} catch (IOException e) {
			System.err.println("Pheno file initialization exception.");
			e.printStackTrace(System.err);
		}
	}

	public PedigreeFile getPedigreeData() {
		return pedData;
	}
	
	public GMDRPhenoFile getPhenotypeData() {
		return phenoData;
	}
	
	public MapFile getMapData() {
		return mapData;
	}
}
