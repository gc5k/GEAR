package gear.qc.snpqc;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;

import gear.ConstValues;
import gear.family.pedigree.file.MapFile;
import gear.family.pedigree.file.SNP;
import gear.subcommands.CommandArguments;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;

public class SNPFilter {
	protected MapFile mapData;
	private int[] workingSNPs;
	private boolean[] snpWorkingFlags;

	public SNPFilter(MapFile mapData) {
		this.mapData = mapData;
	}

	public void filter(CommandArguments cmdArgs) {
		snpWorkingFlags = new boolean[mapData.getMarkerNumber()];
		if (cmdArgs.isExtractFile()) {
			selectSNPs(cmdArgs);
		} else if (cmdArgs.isExcludeFile()) {
			removeSNPs(cmdArgs);
		} else if (cmdArgs.isChr()) {
			selectChromosomes(cmdArgs);
		} else if (cmdArgs.isNotChr()) {
			removeChromosomes(cmdArgs);
		} else {
			Arrays.fill(snpWorkingFlags, true);
		}
		collectWorkingSNPs();
		if (workingSNPs.length == 0) {
			Logger.printUserLog("No SNPs were remained for analysis. GEAR quit.");
			System.exit(1);
		}
		Logger.printUserLog(workingSNPs.length + " SNPs were remained for analysis.");
	}

	private void selectSNPs(CommandArguments cmdArgs) {
		BufferedReader eFile = FileUtil.FileOpen(cmdArgs.getExtractFile());
		HashSet<String> snpChosenSet = NewIt.newHashSet();
		String line;
		try {
			while ((line = eFile.readLine()) != null) {
				String[] s = line.split(ConstValues.WHITESPACE_DELIMITER);
				for(int i = 0; i < s.length; i++) snpChosenSet.add(s[i]);
			}
		} catch (IOException e) {
			Logger.handleException(e, "An exception occurred when reading '" + cmdArgs.getExtractFile() + "'.");
		}

		if (snpChosenSet.size() > 0) {
			Logger.printUserLog("Read " + snpChosenSet.size() + " SNPs from '" + cmdArgs.getExtractFile()+ "'.");
			for (int i = 0; i < mapData.getMarkerList().size(); i++) {
				SNP snp = mapData.getMarkerList().get(i);
				String snpName = snp.getName();
				snpWorkingFlags[i] = snpChosenSet.contains(snpName);
			}
		}
	}

	private void removeSNPs(CommandArguments cmdArgs) {
		BufferedReader eFile = FileUtil.FileOpen(cmdArgs.getExcludeFile());
		HashSet<String> snpExcludeSet = NewIt.newHashSet();
		String line;
		try {
			while ((line = eFile.readLine()) != null) {
				String[] s = line.split(ConstValues.WHITESPACE_DELIMITER);
				for (int i = 0; i < s.length; i++) snpExcludeSet.add(s[i]);
			}
		} catch (IOException e) {
			Logger.handleException(e, "An exception occurred when reading '" + cmdArgs.getExcludeFile() + "'.");
		}

		if (snpExcludeSet.size() > 0) {
			Logger.printUserLog("Read " + snpExcludeSet.size() + " SNPs from '" + cmdArgs.getExcludeFile() + "'.");
			for (int i = 0; i < mapData.getMarkerList().size(); i++) {
				SNP snp = mapData.getMarkerList().get(i);
				String snpName = snp.getName();
				snpWorkingFlags[i] = !snpExcludeSet.contains(snpName);
			}
		}
	}

	private void removeChromosomes(CommandArguments cmdArgs) {
		Logger.printUserLog("Filtering out SNPs not at the selected chromosome(s)...");
		HashSet<String> chrNotSet = cmdArgs.getNotChr();
		for (int i = 0; i < mapData.getMarkerList().size(); i++) {
			SNP snp = mapData.getMarkerList().get(i);
			String chr = snp.getChromosome();
			snpWorkingFlags[i] = !chrNotSet.contains(chr);
		}
	}

	private void selectChromosomes(CommandArguments cmdArgs) {
		HashSet<String> chrSet = cmdArgs.getChr();
		for (int i = 0; i < mapData.getMarkerList().size(); i++) {
			SNP snp = mapData.getMarkerList().get(i);
			String chr = snp.getChromosome();
			snpWorkingFlags[i] = chrSet.contains(chr);
		}
	}

	private void collectWorkingSNPs() {
		ArrayList<Integer> workingSNPsArrayList = new ArrayList<Integer>();
		for (int i = 0; i < snpWorkingFlags.length; ++i) {
			if (snpWorkingFlags[i])
				workingSNPsArrayList.add(i);
		}
		workingSNPs = new int[workingSNPsArrayList.size()];
		for (int i = 0; i < workingSNPs.length; ++i) {
			workingSNPs[i] = workingSNPsArrayList.get(i);
		}
	}
	
	public boolean isSnpIncluded(int snpIndex) {
		return snpWorkingFlags[snpIndex];
	}

	public int[] getWorkingSNPs() {
		return workingSNPs;
	}
}
