package gear.family.qc.colqc;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;

import gear.ConstValues;
import gear.family.pedigree.file.MapFile;
import gear.family.pedigree.file.SNP;
import gear.subcommands.CommandArguments;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;

public class SNPFilter implements SNPFilterInterface {

	protected MapFile mapData;
	protected int[] WSNP;

	protected int[] wseq = null;
	protected int[] bgseq = null;
	private ArrayList<SNP> snpList;
	private HashSet<String> snpChosenSet;

	protected HashSet<Integer> selectedSNPSet;
	protected HashSet<Integer> excludedSNPSet;

	public SNPFilter(MapFile mapData) {
		this.mapData = mapData;
		snpList = mapData.getMarkerList();
		selectedSNPSet = NewIt.newHashSet();
		excludedSNPSet = NewIt.newHashSet();
	}

	public void SelectSNP() {
		makeWSNPList();
		return;
	}

	public void SelectSNP(CommandArguments cmdArgs) {

		if (cmdArgs.isExtractFile()) {
			selectSNP(cmdArgs);
		} else if (cmdArgs.isExcludeFile()) {
			removeSNP(cmdArgs);
		} else if (cmdArgs.isChr()) {
			selectChromosome(cmdArgs);
		} else if (cmdArgs.isNotChr()) {
			removeChromosome(cmdArgs);
		}

		makeWSNPList();
		if (WSNP.length == 0) {
			Logger.printUserLog("No SNPs were remained for analysis. GEAR quit.");
			System.exit(1);
		}
		Logger.printUserLog(WSNP.length + " SNPs were remained for analysis.");
		return;
	}

	private void selectSNP(CommandArguments cmdArgs) {
		BufferedReader eFile = FileUtil.FileOpen(cmdArgs.getExtractFile());
		snpChosenSet = NewIt.newHashSet();
		String line;
		try {
			while ((line = eFile.readLine()) != null) {
				String[] s = line.split(ConstValues.WHITESPACE_DELIMITER);
				for(int i = 0; i < s.length; i++) snpChosenSet.add(s[i]);
			}
		} catch (IOException e) {
			Logger.handleException(e, "An exception occurred when reading '"
					+ cmdArgs.getExtractFile() + "'.");
		}

		if (snpChosenSet.size() > 0) {
			Logger.printUserLog("Read " + snpChosenSet.size() + " SNPs from '"
					+ cmdArgs.getExtractFile()+ "'.");
			for (int i = 0; i < snpList.size(); i++) {
				SNP snp = snpList.get(i);
				String snpName = snp.getName();
				if (snpChosenSet.contains(snpName)) {
					includeSNP(i);
				}
			}
		}
	}

	private void removeSNP(CommandArguments cmdArgs) {

		BufferedReader eFile = FileUtil.FileOpen(cmdArgs.getExcludeFile());
		snpChosenSet = NewIt.newHashSet();
		String line;
		try {
			while ((line = eFile.readLine()) != null) {
				String[] s = line.split(ConstValues.WHITESPACE_DELIMITER);
				for (int i = 0; i < s.length; i++) snpChosenSet.add(s[i]);
			}
		} catch (IOException e) {
			Logger.handleException(e, "An exception occurred when reading '"
					+ cmdArgs.getExcludeFile() + "'.");
		}

		if (snpChosenSet.size() > 0) {
			Logger.printUserLog("Read " + snpChosenSet.size() + " SNPs from '"
					+ cmdArgs.getExcludeFile() + "'.");
			for (int i = 0; i < snpList.size(); i++) {
				SNP snp = snpList.get(i);
				String snpName = snp.getName();
				if (snpChosenSet.contains(snpName)) {
					excludeSNP(i);
				}
			}
		}
	}

	private void removeChromosome(CommandArguments cmdArgs) {

		Logger.printUserLog("Filtering out SNPs not at the selected chromosome(s)...");
		HashSet<String> chrNotSet = cmdArgs.getNotChr();
		for (int i = 0; i < snpList.size(); i++) {
			SNP snp = snpList.get(i);
			String chr = snp.getChromosome();
			if (chrNotSet.contains(chr)) {
				excludeSNP(i);
			} else {
				includeSNP(i);
			}
		}
	}

	private void selectChromosome(CommandArguments cmdArgs) {

		Logger.printUserLog("Choosing SNPs from the selected chromosome(s)...");
		HashSet<String> chrSet = cmdArgs.getChr();
		for (int i = 0; i < snpList.size(); i++) {
			SNP snp = snpList.get(i);
			String chr = snp.getChromosome();
			if (chrSet.contains(chr)) {
				includeSNP(i);
			} else {
				excludeSNP(i);
			}
		}
	}

	private void makeWSNPList() {

		if (selectedSNPSet.size() > 0) {
			WSNP = new int[selectedSNPSet.size()];
			int c = 0;
			for (Iterator<Integer> e = selectedSNPSet.iterator(); e.hasNext();) {
				Integer V = e.next();
				WSNP[c++] = V.intValue();
			}
			Arrays.sort(WSNP);
		} else if (excludedSNPSet.size() > 0) {
			WSNP = new int[snpList.size() - excludedSNPSet.size()];
			int c = 0;
			for (int i = 0; i < snpList.size(); i++) {
				if (!excludedSNPSet.contains(i)) {
					WSNP[c++] = i;
				}
			}
		} else {
			WSNP = new int[snpList.size()];
			for (int i = 0; i < snpList.size(); i++) {
				WSNP[i] = i;
			}
		}

		wseq = new int[WSNP.length];
		for (int i = 0; i < WSNP.length; i++) {
			wseq[i] = i;
		}
	}

	private boolean includeSNP(int i) {
		if (selectedSNPSet.contains(i)) {
			return true;
		} else {
			selectedSNPSet.add(i);
			return false;
		}
	}

	private boolean excludeSNP(int i) {

		if (excludedSNPSet.contains(i)) {
			return true;
		} else {
			excludedSNPSet.add(i);
			return false;
		}
	}

	public int[] getWorkingSNP() {
		return WSNP;
	}

	@Override
	public int[] getWSeq() {
		return wseq;
	}

	@Override
	public int[] getBgSeq() {
		return bgseq;
	}

	@Override
	public int[][] getWSeq2() {
		return null;
	}
}
