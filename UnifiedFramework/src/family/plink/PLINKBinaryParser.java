package family.plink;

import java.io.IOException;

import test.Test;

import family.pedigree.file.BEDReader;
import family.pedigree.file.BIMReader;
import family.pedigree.file.PhenotypeFile;

public class PLINKBinaryParser extends PLINKParser {

	protected String FamFile;
	public PLINKBinaryParser(String ped, String map, String fam, String phe) {
		super(ped, map, phe);
		FamFile = fam;
	}

	@Override
	public void Parse() {
		mapData = new BIMReader(mapFile);

		if (phenotypeFile != null) {
			phenoData = new PhenotypeFile();
			ParsePhenoFile();
		}
		if (mapFile != null) {
			ParseMapFile();
			Test.LOG.append(mapData.getMarkerNumber() + " markers.");
			Test.LOG.append("\n");
			System.err.println(mapData.getMarkerNumber() + " markers.");
			pedData = new BEDReader(FamFile, snpFilter.getWorkingSNP().length, mapData);
			pedData.setHeader(false);
			ParsePedFile();
			Test.LOG.append(pedData.getNumIndividuals() + " individuals.");
			Test.LOG.append("\n");
			System.err.println(pedData.getNumIndividuals() + " individuals.");
		}
//		mapData.setPolymorphism(pedData.getAlleleFrequency());
		pedData.cleanup();
		if (phenoData != null) {
			Test.LOG.append(phenoData.getNumTraits() + " traits.");
			Test.LOG.append("\n");
			System.err.println(phenoData.getNumTraits() + " traits.");
			
		}
	}

	@Override
	public void ParsePedFile() {
		try {
			pedData.parseLinkage(pedigreeFile, mapData.getMarkerNumberOriginal(), snpFilter.getWorkingSNP());
		} catch (IOException e) {
			System.err.println("Pedgree file initialization exception.");
			e.printStackTrace(System.err);
		}
	}
}
