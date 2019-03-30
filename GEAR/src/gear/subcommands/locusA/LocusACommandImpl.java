package gear.subcommands.locusA;

import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;

import gear.family.GenoMatrix.GenotypeMatrix;
import gear.family.pedigree.Hukou;
import gear.family.pedigree.file.BEDReader;
import gear.family.pedigree.file.MapFile;
import gear.family.pedigree.file.SNP;
import gear.family.plink.PLINKBinaryParser;
import gear.family.plink.PLINKParser;
import gear.qc.sampleqc.SampleFilter;
import gear.qc.snpqc.SNPFilter;
import gear.qc.snpqc.SNPFilterPostQC;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.subcommands.locus.LocusCommandArguments;
import gear.util.FileUtil;
import gear.util.Logger;

public class LocusACommandImpl extends CommandImpl {
	private LocusACommandArguments locusArgs;
	private MapFile map;
	private PrintStream resultFile;
	private DecimalFormat fmt1 = new DecimalFormat("0.0000");
	private DecimalFormat fmt2 = new DecimalFormat("0.00E000");

	private SampleFilter sf;
	private GenotypeMatrix pGM;

	@Override
	public void execute(CommandArguments cmdArgs) {
		
		locusArgs = (LocusACommandArguments) cmdArgs;

		PLINKParser pp = PLINKParser.parse(locusArgs);
		sf = new SampleFilter(pp.getPedigreeData(), cmdArgs);

		pGM = new GenotypeMatrix(sf.getSample(), pp.getMapData(), cmdArgs);

		resultFile = FileUtil.CreatePrintStream(locusArgs.getOutRoot() + ".locus");
		resultFile.println("SNP\tCHR\tBP\tRefAllele\tAltAllele\tFreq\tVar\tEVar\tAA\tAa\taa\tnChr");

		executeBedSnpMajor();
		return;
	}
	
	private void executeBedSnpMajor() {
		float[][] locusStat = pGM.getQCedAlleleStat();
		for (int i = 0; i < pGM.getNumMarker(); ++i) {
			SNP snp = pGM.getSNPList().get(i);
			printResultOfSNP(snp, locusStat[i][0], locusStat[i][2], 
					locusArgs.isInbred()? 4* locusStat[i][0]*locusStat[i][1]: 2* locusStat[i][0]*locusStat[i][1], 
							locusStat[i][3], locusStat[i][4], locusStat[i][5]);
		}
		resultFile.close();
		Logger.printUserLog("Saved " + pGM.getNumMarker() + " results to " + locusArgs.getOutRoot() + ".locus.");
	}
	
	private void printResultOfSNP(
			SNP snp,
			float alleleFreq0,
			float alleleVar,
			float eVar,
			float genoCnt_AA,
			float genoCnt_Aa,
			float genoCnt_aa) {
		if (Float.isNaN(alleleFreq0)) {
			resultFile.println(snp.getName() + "\t" + snp.getChromosome() + "\t" + snp.getPosition() + "\t"
					+ snp.getFirstAllele() + "\t" + snp.getSecAllele() + "\t"
					+ "NA" + "\t"
					+ "NA" + "\t"
					+ "NA" + "\t"
					+ (int) genoCnt_AA + "\t"
					+ (int) genoCnt_Aa + "\t"
					+ (int) genoCnt_aa + "\t"
					+ (((int) genoCnt_AA + (int) genoCnt_Aa + (int) genoCnt_aa) << 1));			
			
		} else {
			resultFile.println(snp.getName() + "\t" + snp.getChromosome() + "\t" + snp.getPosition() + "\t"
					+ snp.getFirstAllele() + "\t" + snp.getSecAllele() + "\t"
					+ (alleleFreq0 > 0.0001 ? fmt1.format(alleleFreq0) : fmt2.format(alleleFreq0)) + "\t"
					+ (alleleVar > 0.0001 ? fmt1.format(alleleVar) : fmt2.format(alleleVar)) + "\t"
					+ (eVar > 0.001 ? fmt1.format(eVar) : fmt2.format(eVar)) + "\t"
					+ (int) genoCnt_AA + "\t"
					+ (int) genoCnt_Aa + "\t"
					+ (int) genoCnt_aa + "\t"
					+ (((int) genoCnt_AA + (int) genoCnt_Aa + (int) genoCnt_aa) << 1));			
		}
	}
}
