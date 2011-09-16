package family.plink;

import java.io.IOException;

import admixture.parameter.Parameter;
import family.pedigree.file.MapFile;
import family.pedigree.file.PhenotypeFile;
import family.pedigree.file.TransposePedigreeReader;

public class PLINKTransposeParser extends PLINKParser {

	protected String FamFile;
	public PLINKTransposeParser(String ped, String Fam, String map, String phe) {
		super(ped, null, phe);
		FamFile = Fam;
		// TODO Auto-generated constructor stub
	}

	public void Parse() {
		mapData = new MapFile(null);
		pedData = new TransposePedigreeReader(pedigreeFile, FamFile, mapData);
		pedData.setHeader(Parameter.header);
		if (phenotypeFile != null) {
			phenoData = new PhenotypeFile();
			ParsePhenoFile();
		}
		pedData.setHeader(false);
		ParsePedFile();

//		mapData.setPolymorphism(pedData.getPolymorphism(), pedData.getAlleleFrequency());
		mapData.setPolymorphismMarker(pedData.getPolymorphism());
		pedData.cleanup();
	}

	public void ParsePedFile() {

		try {
			pedData.parseLinkage(pedigreeFile, 0, snpFilter.getWorkingSNP());
		} catch (IOException e) {
			System.err.println("Pedgree file initialization exception.");
			e.printStackTrace(System.err);
		}
	}
}
