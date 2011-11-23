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
			Test.LOG.append("reading " + phenotypeFile + ".\n");
			System.err.println("reading " + phenotypeFile + ".");
			Test.LOG.append(phenoData.getNumTraits() + " traits in " + phenotypeFile + ".\n");
			System.err.println(phenoData.getNumTraits() + " traits in " + phenotypeFile + ".");

		}
		if (mapFile != null) {
			ParseMapFile();
			Test.LOG.append("reading " + mapFile + ".\n");
			Test.LOG.append(mapData.getMarkerNumberOriginal() + " markers in " + mapFile + ".\n");
//			Test.LOG.append(mapData.getMarkerNumber() + " selected markers.\n");
			System.err.println("reading " + mapFile + ".");
			System.err.println(mapData.getMarkerNumberOriginal() + " markers in " + mapFile + ".");
//			System.err.println(mapData.getMarkerNumber() + " selected markers.");
			pedData = new BEDReader(FamFile, snpFilter.getWorkingSNP().length, mapData);
			pedData.setHeader(false);
			ParsePedFile();
			Test.LOG.append("reading " + pedigreeFile + ".");
			Test.LOG.append(pedData.getNumIndividuals() + " individuals.\n");
			System.err.println("reading " + pedigreeFile + ".");
			System.err.println(pedData.getNumIndividuals() + " individuals.");
		}
//		mapData.setPolymorphism(pedData.getAlleleFrequency());
		pedData.cleanup();
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
