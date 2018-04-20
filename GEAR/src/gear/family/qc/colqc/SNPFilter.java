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

public class SNPFilter implements SNPFilterInterface
{

	protected MapFile mapData;
	protected int[] WSNP;

	protected int[] wseq = null;
	protected int[] bgseq = null;
	ArrayList<SNP> snpList;
	HashSet<String> snpChosenSet;

	protected HashSet<Integer> selectedSNPSet;
	protected HashSet<Integer> excludedSNPSet;
	protected HashSet<Integer> bgSNPSet;

	protected int[][] in_range = null;
	protected int[][] ex_range = null;

	public SNPFilter(MapFile mapData)
	{
		this.mapData = mapData;
		snpList = mapData.getMarkerList();
		selectedSNPSet = NewIt.newHashSet();
		excludedSNPSet = NewIt.newHashSet();
		bgSNPSet = NewIt.newHashSet();
	}

	public void Select()
	{
		makeWSNPList();
		return;
	}

	public void Select(CommandArguments cmdArgs)
	{
		if (cmdArgs.isChr()  || cmdArgs.isNotChr()) {
			selectChromosome(cmdArgs);
		}

		if (cmdArgs.isExtractFile() || cmdArgs.isExcludeFile()) {
			chooseSNP(cmdArgs);
		}

		makeWSNPList();
		Logger.printUserLog(WSNP.length + " SNPs were included for analysis.");			
		Logger.printUserLog(excludedSNPSet.size() + " SNPs were removed from analysis.");			
		return;
	}

	private void chooseSNP(CommandArguments cmdArgs)
	{
		boolean isExtract = cmdArgs.getExtractFile() != null ? true:false;
		BufferedReader eFile = FileUtil.FileOpen(isExtract? cmdArgs.getExtractFile():cmdArgs.getExcludeFile());
		snpChosenSet = NewIt.newHashSet();
		String line;
		try {
			while ((line = eFile.readLine()) != null) {
				String[] s = line.split(ConstValues.WHITESPACE_DELIMITER);
				snpChosenSet.add(s[0]);
			}
		} 
		catch (IOException e) {
			Logger.handleException(e,
					"An exception occurred when reading '"
							+ (isExtract?cmdArgs.getExtractFile():cmdArgs.getExcludeFile())
							+ "'.");
		}

		if (snpChosenSet.size() > 0) {
			Logger.printUserLog("Read " + snpChosenSet.size() + " SNPs from '" + (isExtract?cmdArgs.getExtractFile():cmdArgs.getExcludeFile()) + "'.");
			for (int i = 0; i < snpList.size(); i++) {
				SNP snp = snpList.get(i);
				String snpName = snp.getName();
				if (snpChosenSet.contains(snpName)) {
					if (isExtract) {
						includeSNP(i);
					} else {
						excludeSNP(i);
					}
				}
			}
		}
	}

	private void selectChromosome(CommandArguments cmdArgs)
	{

		if (cmdArgs.isChr()) {
			Logger.printUserLog("Matching SNPs to the selected chromosome(s)...");
			HashSet<String> chrSet = cmdArgs.getChr();
			for (int i = 0; i < snpList.size(); i++)	{
				SNP snp = snpList.get(i);
				String chr = snp.getChromosome();
				if (chrSet.contains(chr)) {
					includeSNP(i);
				} else {
					excludeSNP(i);
				}
			}
		}

		if (cmdArgs.isNotChr()) {
			Logger.printUserLog("Filtering out SNPs at not selected chromosome(s)...");
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
	}

	private void makeWSNPList()
	{

		if (selectedSNPSet.size() > 0)
		{
			WSNP = new int[selectedSNPSet.size()];
			int c = 0;
			for (Iterator<Integer> e = selectedSNPSet.iterator(); e.hasNext();)
			{
				Integer V = e.next();
				WSNP[c++] = V.intValue();
			}
			Arrays.sort(WSNP);
		} else if (excludedSNPSet.size() > 0)
		{
			WSNP = new int[snpList.size() - excludedSNPSet.size()];
			int c = 0;
			for (int i = 0; i < snpList.size(); i++)
			{
				if (!excludedSNPSet.contains(i))
				{
					WSNP[c++] = i;
				}
			}
		} else
		{
			WSNP = new int[snpList.size()];
			for (int i = 0; i < snpList.size(); i++)
			{
				WSNP[i] = i;
			}
		}

		wseq = new int[WSNP.length];
		for (int i = 0; i < WSNP.length; i++)
		{
			wseq[i] = i;
		}
	}

	private boolean includeSNP(int i)
	{
		if (selectedSNPSet.contains(i))
		{
			return true;
		} else
		{
			selectedSNPSet.add(i);
			return false;
		}
	}

	private boolean excludeSNP(int i)
	{

		if (excludedSNPSet.contains(i))
		{
			return true;
		} else
		{
			excludedSNPSet.add(i);
			return false;
		}
	}

	public int[] getWorkingSNP()
	{
		return WSNP;
	}

	@Override
	public int[] getWSeq()
	{
		return wseq;
	}

	@Override
	public int[] getBgSeq()
	{
		return bgseq;
	}

	@Override
	public int[][] getWSeq2()
	{
		return null;
	}
}
