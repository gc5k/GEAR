package family.plink;

import java.io.File;
import java.io.IOException;

import family.pedigree.file.GMDRPhenoFile;
import family.pedigree.file.GMDRPhenoFileException;
import family.pedigree.file.MDRPedFileException;
import family.pedigree.file.MapFile;
import family.pedigree.file.PedigreeFile;

public class PLINKParser {
	protected MapFile MapData;
	protected PedigreeFile PedData;
	protected GMDRPhenoFile PhenoData;

	protected String pedigreeFile;
	protected String phenotypeFile;
	protected String mapFile;
	
	public PLINKParser (String ped, String phe, String map) {
		pedigreeFile = ped;
		phenotypeFile = phe;
		mapFile = map;

		MapData = new MapFile(mapFile);
		PhenoData = new GMDRPhenoFile();
		PedData = new PedigreeFile();
		initial();
	}

	private void initial() {
		if(mapFile != null) {
			ParseMapFile();
			PedData.setHeader(false);
			ParsePedFile();
		} else {
			PedData.setHeader(true);
			ParsePedFile();
			MapData.setMarker(PedData.getNumMarker());
		}

		if (phenotypeFile != null) {
			ParsePhenoFile();
		}
	}

	public void ParseMapFile() {
		if(mapFile != null) {
			MapData.parseMap();
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
			PedData.Initial(PedFile);
			PedData.parseLinkage();
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
			PhenoData.Initial(PheFile);
			PhenoData.parsePhenotype();
		} catch (GMDRPhenoFileException e) {
			System.err.println("Pheno file initialization exception.");
			e.printStackTrace(System.err);
		} catch (IOException e) {
			System.err.println("Pheno file initialization exception.");
			e.printStackTrace(System.err);
		}
	}

}
