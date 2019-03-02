package gear.subcommands.locus;

import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;

import gear.family.pedigree.Hukou;
import gear.family.pedigree.file.BEDReader;
import gear.family.pedigree.file.MapFile;
import gear.family.pedigree.file.SNP;
import gear.family.plink.PLINKBinaryParser;
import gear.family.plink.PLINKParser;
import gear.qc.sampleqc.SampleFilter;
import gear.qc.snpqc.SNPFilter;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.FileUtil;
import gear.util.Logger;

public class LocusCommandImpl extends CommandImpl {
	private LocusCommandArguments locusArgs;
	private MapFile map;
	private SNPFilter snpFilter;
	private BEDReader bed;
	private PrintStream resultFile;
	private DecimalFormat fmt1 = new DecimalFormat("0.0000");
	private DecimalFormat fmt2 = new DecimalFormat("0.00E000");

	@Override
	public void execute(CommandArguments cmdArgs) {
		locusArgs = (LocusCommandArguments) cmdArgs;

		PLINKParser pp = PLINKParser.create(locusArgs);
		pp.parseSmallFiles();
		map = pp.getMapData();
		snpFilter = pp.getSNPFilter();

		SampleFilter samFilter = new SampleFilter(pp.getPedigreeData(), locusArgs);

		resultFile = FileUtil.CreatePrintStream(locusArgs.getOutRoot() + ".locus");
		resultFile.println("SNP\tCHR\tBP\tRefAllele\tAltAllele\tFreq\tVar\tEVar\tAA\tAa\taa\tnChr");

		bed = (BEDReader) pp.getPedigreeData();
		executeBedSnpMajor();
		return;

	}
	
	private void executeBedSnpMajor() {
		int numMarkers = map.getMarkerNumberOriginal();
		int numSamples = bed.getNumIndividuals();
		int workingSnpIndex = 0;
		
		ArrayList<Hukou> hkBook = bed.getHukouBook();

		for (int i = 0; i < numMarkers; ++i) {
			if (snpFilter.isSnpIncluded(i)) {
				int genoCnt_AA = 0;
				int genoCnt_Aa = 0;
				int genoCnt_aa = 0;
				int missingCnt = 0;
				int sum = 0;
				int squareSum = 0;
				for (int j = 0; j < numSamples; j += 4) {
					int nextByte = bed.readNextByte();
					int indCnt = j;
					for (int k = 0; k < 8; k += 2) {

						if(indCnt == numSamples) break;

						if (!hkBook.get(indCnt++).isAvailable()) {
							continue;
						}

						int genotype = (nextByte >> k) & 0b11;

						switch (genotype) {
						case PLINKBinaryParser.HOMOZYGOTE_FIRST:
							++genoCnt_AA;
							break;
						case PLINKBinaryParser.HETEROZYGOTE:
							++genoCnt_Aa;
							sum += 1;
							squareSum += 1;
							break;
						case PLINKBinaryParser.HOMOZYGOTE_SECOND:
							++genoCnt_aa;
							sum += 2;
							squareSum += 4;
							break;
						case PLINKBinaryParser.MISSING_GENOTYPE:
							++missingCnt;
							break;
						}
					}
				}
				int validSampleCnt = numSamples - missingCnt;
				double variance = 0;
				double alleleFreq0 = 0;
				double alleleFreq1 = 0;

				if (validSampleCnt > 2) {
					alleleFreq1 = (double)sum / (2.0*validSampleCnt);
					alleleFreq0 = 1-alleleFreq1;
					variance = (squareSum - validSampleCnt * alleleFreq1 * alleleFreq1) / (validSampleCnt - 1);
				}
				double eVariance = calculateEVariance(alleleFreq0, alleleFreq1);

				if (locusArgs.isMAF() || locusArgs.isMaxMAF()
						|| locusArgs.isGENO() || locusArgs.isMAFRange()) {
					double maf = alleleFreq0<alleleFreq1 ? alleleFreq0 : alleleFreq1;
					if (locusArgs.isMAF() && maf < locusArgs.getMAF()) {
						continue;
					}
					if (locusArgs.isMaxMAF() && maf > locusArgs.getMaxMAF()) {
						continue;
					}
					if (locusArgs.isGENO() && (missingCnt * 1.0D/numSamples) > locusArgs.getGENO()) {
						continue;
					}
					if (locusArgs.isMAFRange() ) {
						double[][] mafR = locusArgs.getMAFRange();
						boolean mafFlag = false;
						for (int j = 0; j < mafR.length; j++) {
							if (maf >= mafR[j][0] && maf <= mafR[j][1]) {
								mafFlag = true;
							}
						}
						if (!mafFlag) {
							continue;
						}
					}
				}
				SNP snp = map.getSNP(workingSnpIndex++);
				printResultOfSNP(snp, alleleFreq0, variance, eVariance, genoCnt_AA, genoCnt_Aa, genoCnt_aa);
			} else {
				bed.skipOneRow();
			}
		}
		resultFile.close();
		Logger.printUserLog("Save " + workingSnpIndex + " results to " + locusArgs.getOutRoot() + ".locus.");
	}
	
	private double calculateEVariance(double alleleFreq0, double alleleFreq1) {
		return locusArgs.isInbred() ? 4 * alleleFreq0 * alleleFreq1 : 2 * alleleFreq0 * alleleFreq1;
	}
	
	private void printResultOfSNP(
			SNP snp,
			double alleleFreq0,
			double alleleVar,
			double eVar,
			int genoCnt_AA,
			int genoCnt_Aa,
			int genoCnt_aa) {
		resultFile.println(snp.getName() + "\t" + snp.getChromosome() + "\t" + snp.getPosition() + "\t"
				+ snp.getFirstAllele() + "\t" + snp.getSecAllele() + "\t"
				+ (alleleFreq0 > 0.0001 ? fmt1.format(alleleFreq0) : fmt2.format(alleleFreq0)) + "\t"
				+ (alleleVar > 0.0001 ? fmt1.format(alleleVar) : fmt2.format(alleleVar)) + "\t"
				+ (eVar > 0.001 ? fmt1.format(eVar) : fmt2.format(eVar)) + "\t"
				+ genoCnt_AA + "\t"
				+ genoCnt_Aa + "\t"
				+ genoCnt_aa + "\t"
				+ ((genoCnt_AA + genoCnt_Aa + genoCnt_aa) << 1));
	}
}
