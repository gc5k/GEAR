package family.mdr.filter;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;

import org.apache.commons.lang3.ArrayUtils;

import admixture.parameter.Parameter;
import family.pedigree.file.MapFile;
import family.pedigree.file.SNP;

import test.Test;
import util.NewIt;

public class SNPFilterI implements SNPFilterInterface {
	protected MapFile mapData;
	protected int[] WSNP;

	protected int[][] wseq = null;
	protected int[] bgseq = null;
	protected boolean filterFlag = false;
	ArrayList<SNP> snpList;

	protected HashSet<Integer> selectedSNPSet;
	protected HashSet<Integer> excludedSNPSet;
	protected HashSet<Integer> bgSNPSet;

	protected ArrayList<HashSet<Integer>> snpArrays;

	protected int[][] in_range = null;
	protected int[][] ex_range = null;

	protected int[][] chrSNP = null;

	protected int[][] snpWin = null;

	public SNPFilterI(MapFile mapData) {
		this.mapData = mapData;
		snpList = mapData.getMarkerList();
		selectedSNPSet = NewIt.newHashSet();
		excludedSNPSet = NewIt.newHashSet();
		bgSNPSet = NewIt.newHashSet();
		snpArrays = NewIt.newArrayList();
	}

	public void Select() {

		if (Parameter.bgsnpFlag) {
			selectBackgroundSNP();
			filterFlag = true;
		}

		if (Parameter.inchrFlag) {
			selectChromosome();
			filterFlag = true;
		}

		if (Parameter.snpwindowFlag) {
			selectSNPWindow();
			filterFlag = true;
		}

		if (Parameter.snpFlag) {
			selectSNPs();
			filterFlag = true;
		}

		if (Parameter.regionFlag) {
			selectSNPRegion();
			filterFlag = true;
		}

		if (Parameter.geneFlag) {
			selectGene();
			filterFlag = true;
		}

		if (filterFlag && selectedSNPSet.size() == 0
				&& excludedSNPSet.size() == 0 && bgSNPSet.size() == 0) {
			System.err.println("No snps selected. GMDR quit.");
			Test.LOG.append("No snps selected. GMDR quit.\n");
			Test.printLog();
			System.exit(0);
		}
		makeWSNPList();

		return;
	}

	private void selectBackgroundSNP() {

		for (int i = 0; i < snpList.size(); i++) {
			SNP snp = snpList.get(i);
			String chr = snp.getName();
			int idx = ArrayUtils.indexOf(Parameter.bgsnp, chr);
			if (idx >= 0) {
				addBackgroundSNP(i);
			}
		}

	}

	private void selectChromosome() {

		ArrayList<HashSet<Integer>> chrSNPSet = NewIt.newArrayList();
		int L = Parameter.in_chr.length;
		for (int i = 0; i < L; i++) {
			HashSet<Integer> chrSet = NewIt.newHashSet();
			chrSNPSet.add(chrSet);
		}

		for (int i = 0; i < snpList.size(); i++) {
			SNP snp = snpList.get(i);
			String chr = snp.getChromosome();
			int idx = ArrayUtils.indexOf(Parameter.in_chr, chr);
			if (idx >= 0) {
				includeSNP(i);
				HashSet<Integer> chrSet = chrSNPSet.get(idx);
				chrSet.add(new Integer(i));
			}
		}

		chrSNP = new int[L][];
		for (int i = 0; i < L; i++) {
			HashSet<Integer> chrSet = chrSNPSet.get(i);
			if (chrSet.size() == 0) {
				System.err.println("could not find chromosome "
						+ Parameter.in_chr[i]);
				continue;
			}
			snpArrays.add(chrSet);
			chrSNP[i] = new int[chrSet.size()];
			int c = 0;
			for (Iterator<Integer> e = chrSet.iterator(); e.hasNext();) {
				chrSNP[i][c++] = e.next().intValue();
			}
		}
	}

	private void selectSNPWindow() {
		ArrayList<SNP> snpwindowList = NewIt.newArrayList();
		ArrayList<Integer> snpwindowIndex = NewIt.newArrayList();
		ArrayList<HashSet<Integer>> snpwindowSet = NewIt.newArrayList();

		for (int i = 0; i < Parameter.snpwindow.length; i++) {
			HashSet<Integer> sw = NewIt.newHashSet();
			snpwindowSet.add(sw);
		}

		int c = 0;
		for (int i = 0; i < snpList.size(); i++) {
			SNP snp = snpList.get(i);
			String rs = snp.getName();
			int idx = ArrayUtils.indexOf(Parameter.snpwindow, rs);
			if (idx >= 0) {
				snpwindowList.add(snp);
				snpwindowIndex.add(new Integer(idx));
				c++;
			}
			if (c == Parameter.snpwindow.length) {
				break;
			}
		}
		if (snpwindowList.size() == 0) {
			System.err.println("could not find the snp windows specified.");

		}
		for (int i = 0; i < snpList.size(); i++) {
			SNP snp = snpList.get(i);
			for (int j = 0; j < snpwindowList.size(); j++) {
				int d = snpwindowList.get(j).compareTo(snp);
				if (d != Integer.MAX_VALUE && d != Integer.MIN_VALUE) {
					int idx = snpwindowIndex.get(j).intValue();
					if (d > Parameter.snp_window[idx][0]
							&& d < Parameter.snp_window[idx][1]) {
						includeSNP(i);
						snpwindowSet.get(idx).add(new Integer(i));
					}
				}
			}
		}

		snpWin = new int[Parameter.snpwindow.length][];
		for (int i = 0; i < snpwindowSet.size(); i++) {
			HashSet<Integer> swSet = snpwindowSet.get(i);
			if (swSet.size() == 0) {
				System.err
						.println("could not find snps near "
								+ Parameter.snpwindow[snpwindowIndex.get(i)
										.intValue()]);
				continue;
			}
			snpArrays.add(swSet);
			snpWin[i] = new int[snpwindowSet.get(i).size()];
			c = 0;
			for (Iterator<Integer> e = swSet.iterator(); e.hasNext();) {
				snpWin[i][c++] = e.next().intValue();
			}
		}
	}

	private void selectGene() {

		int count = 0;
		ArrayList<ArrayList<Integer>> xsnps = NewIt.newArrayList();

		if (!Parameter.snp2genefileFlag && !Parameter.snp2genefilesFlag) {
			for (int i = 0; i < Parameter.gene.length; i++) {
				ArrayList<Integer> s = NewIt.newArrayList();
				xsnps.add(s);
			}

			for (int i = 0; i < snpList.size(); i++) {
				SNP snp = snpList.get(i);
				String chr = snp.getChromosome();
				int pos = snp.getPosition();
				for (int j = 0; j < Parameter.gene_chr.length; j++) {

					if (chr.compareTo(Parameter.gene_chr[j]) == 0) {
						if (pos >= (Parameter.gene_begin[j] - Parameter.genewindow) * 1000
								&& pos <= (Parameter.gene_end[j] + Parameter.genewindow) * 1000) {
							xsnps.get(j).add(i);
							includeSNP(i);
							count++;
						}
					}
				}
			}
		} else {
			for (int i = 0; i < Parameter.gene.length; i++) {
				ArrayList<Integer> s = NewIt.newArrayList();
				xsnps.add(s);
			}

			for (Iterator<Integer> e = selectedSNPSet.iterator(); e.hasNext();) {
				int i = e.next().intValue();
				SNP snp = snpList.get(i);
				String chr = snp.getChromosome();
				int pos = snp.getPosition();
				for (int j = 0; j < Parameter.gene_chr.length; j++) {
					if (chr.compareTo(Parameter.gene_chr[j]) == 0) {
						if (pos >= (Parameter.gene_begin[j] - Parameter.genewindow) * 1000
								&& pos <= (Parameter.gene_end[j] + Parameter.genewindow) * 1000) {
							xsnps.get(j).add(i);
							count++;
						}
					}
				}
			}
		}

		int countgene = 0;
		if (Parameter.snp2genefilesFlag) {
			for (int i = 0; i < xsnps.size(); i++) {
				StringBuffer sbsnp = new StringBuffer(Parameter.gene[i]
						+ ".snp");
				StringBuffer sb = new StringBuffer(Parameter.gene[i] + ".gene");
				PrintStream PW = null;
				PrintStream PWsnp = null;
				try {
					PW = new PrintStream(sb.toString());
					PWsnp = new PrintStream(sbsnp.toString());
				} catch (FileNotFoundException e) {
					e.printStackTrace();
				}

				ArrayList<Integer> s = xsnps.get(i);
				if (s.size() == 0) {
					continue;
				}
				countgene++;
				System.err.print(s.size() + " snps selected with gene "
						+ Parameter.gene[i]);
				System.err.println(" [writing snps into " + Parameter.gene[i]
						+ ".gene]");
				Test.LOG.append(s.size() + " snps selected with gene "
						+ Parameter.gene[i]);
				Test.LOG.append(" [writing snps into " + Parameter.gene[i]
						+ ".gene]\n");
				for (int j = 0; j < s.size(); j++) {
					SNP snp = snpList.get(s.get(j));
					PWsnp.println(snp.getName());
					PW.println(snp.getName() + " " + snp.getChromosome() + " "
							+ snp.getPosition() + " " + Parameter.gene[i]);
				}
				PWsnp.close();
				PW.close();
			}
			System.err.println(count + " snps selected for " + countgene
					+ " selected genes in total.");
			Test.LOG.append(count + " snps selected for " + countgene
					+ " selected genes in total.\n");
			Test.printLog();
			System.exit(1);
		} else {
			StringBuffer sbsnp = new StringBuffer(Parameter.out + ".snp");
			StringBuffer sb = new StringBuffer(Parameter.out + ".gene");
			PrintStream PW = null;
			PrintStream PWsnp = null;
			try {
				PW = new PrintStream(sb.toString());
				PWsnp = new PrintStream(sbsnp.toString());
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
			for (int i = 0; i < xsnps.size(); i++) {
				ArrayList<Integer> s = xsnps.get(i);
				if (s.size() == 0) {
					continue;
				}
				countgene++;
				System.err.println(s.size() + " snps selected with gene "
						+ Parameter.gene[i]);

				Test.LOG.append(s.size() + " snps selected with gene "
						+ Parameter.gene[i] + "\n");
				for (int j = 0; j < s.size(); j++) {
					SNP snp = snpList.get(s.get(j));
					PWsnp.println(snp.getName());
					PW.println(snp.getName() + " " + snp.getChromosome() + " "
							+ snp.getPosition() + " " + Parameter.gene[i]);
				}
			}
			System.err.println("writing snps into " + Parameter.out + ".gene");
			Test.LOG.append("writing snps into " + Parameter.out + ".gene\n");
			PW.close();
			PWsnp.close();
			System.err.println(count + " snps selected for " + countgene
					+ " selected genes in total.");
			Test.LOG.append(count + " snps selected for " + countgene
					+ " selected genes in total.\n");
			if (Parameter.snp2genefileFlag) {
				Test.printLog();
				System.exit(1);
			}
		}

		for (int i = 0; i < xsnps.size(); i++) {
			ArrayList<Integer> xsnp = xsnps.get(i);
			if (xsnp.size() == 0) {
				continue;
			}
			HashSet<Integer> snpSet = NewIt.newHashSet();
			snpSet.addAll(xsnp);
			snpArrays.add(snpSet);
		}
	}

	private void selectSNPRegion() {
		ArrayList<ArrayList<Integer>> xsnps = NewIt.newArrayList();
		for (int i = 0; i < Parameter.gene.length; i++) {
			ArrayList<Integer> s = NewIt.newArrayList();
			xsnps.add(s);
		}

		for (int i = 0; i < snpList.size(); i++) {
			SNP snp = snpList.get(i);
			String chr = snp.getChromosome();
			int pos = snp.getPosition();
			int idx = ArrayUtils.indexOf(Parameter.chr_reg, chr);
			if (idx >= 0) {
				if (pos >= Parameter.begin[idx] * 1000
						&& pos <= Parameter.end[idx] * 1000) {
					xsnps.get(idx).add(i);
					includeSNP(i);
				}
			}
		}
		for (int i = 0; i < xsnps.size(); i++) {
			ArrayList<Integer> xsnp = xsnps.get(i);
			if (xsnp.size() == 0) {
				continue;
			}
			HashSet<Integer> snpSet = NewIt.newHashSet();
			snpSet.addAll(xsnp);
			snpArrays.add(snpSet);
		}
		return;
	}

	private void selectSNPs() {

		ArrayList<ArrayList<Integer>> xsnps = NewIt.newArrayList();
		for (int i = 0; i < Parameter.xincludesnp.length; i++) {
			ArrayList<Integer> s = NewIt.newArrayList();
			xsnps.add(s);
		}
		for (int i = 0; i < Parameter.xincludesnp.length; i++) {
			HashSet<String> SS = NewIt.newHashSet();

			ArrayList<Integer> s = xsnps.get(i);

			for (int j = 0; j < Parameter.xincludesnp[i].length; j++) {
				SS.add(Parameter.xincludesnp[i][j]);
			}
			for (int j = 0; j < snpList.size(); j++) {
				SNP snp = snpList.get(j);
				String rs = snp.getName();
				if(SS.contains(rs)) {
					includeSNP(j);
					s.add(j);
				}
			}
			System.err.println(s.size() + " of " + SS.size() + " SNPs were selected " + "in " + Parameter.includesnpFile[i] + ".");
			Test.LOG.append(s.size() + " of " + SS.size() + " SNPs were selected " + "in " + Parameter.includesnpFile[i] + ".\n");
		}

		for (int i = 0; i < xsnps.size(); i++) {
			ArrayList<Integer> xsnp = xsnps.get(i);
			if (xsnp.size() == 0) {
				continue;
			}
			HashSet<Integer> snpSet = NewIt.newHashSet();
			snpSet.addAll(xsnp);
			snpArrays.add(snpSet);
		}
		return;
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
		} else {
			if (Parameter.transFlag) {
				System.err.println("no snps selected.");
				Test.LOG.append("no snps selected.\n");
				Test.printLog();
				System.exit(0);
			}
		}

		if (bgSNPSet.size() > 0) {
			bgseq = new int[bgSNPSet.size()];
			int c = 0;
			for (Iterator<Integer> e = bgSNPSet.iterator(); e.hasNext();) {
				int v = e.next().intValue();
				int idx = ArrayUtils.indexOf(WSNP, v);
				bgseq[c++] = idx;
			}
			for (Iterator<Integer> e = bgSNPSet.iterator(); e.hasNext();) {
				Integer V = e.next();
				for (HashSet<Integer> f : snpArrays) {
					if (f.contains(V)) {
						f.remove(V);
					}
				}
			}
		}

		int L = 0;
		for (HashSet<Integer> e : snpArrays) {
			if (e.size() > 0)
				L++;
		}
		if (L == 0) {
			System.err.println("no snps selected.");
			Test.LOG.append("no snps selected.\n");
			Test.printLog();
			System.exit(0);
		}
		wseq = new int[L][];
		int c = 0;
		for (HashSet<Integer> e : snpArrays) {
			if (e.size() > 0) {
				wseq[c] = new int[e.size()];
				int cc = 0;
				for (Integer I : e) {
					int v = I.intValue();
					wseq[c][cc++] = ArrayUtils.indexOf(WSNP, v);
				}
				Arrays.sort(wseq[c]);
				c++;
			}
		}
	}

	private boolean addBackgroundSNP(int i) {
		Integer I = new Integer(i);
		if (bgSNPSet.contains(I)) {
			return true;
		} else {
			bgSNPSet.add(I);
			includeSNP(i);
			return false;
		}
	}

	private boolean includeSNP(int i) {
		Integer I = new Integer(i);
		if (selectedSNPSet.contains(I)) {
			return true;
		} else {
			selectedSNPSet.add(I);
			return false;
		}
	}

	@Override
	public int[] getWorkingSNP() {
		return WSNP;
	}

	@Override
	public int[] getWSeq() {
		return null;
	}

	@Override
	public int[][] getWSeq2() {
		return wseq;
	}

	@Override
	public int[] getBgSeq() {
		return bgseq;
	}
}
