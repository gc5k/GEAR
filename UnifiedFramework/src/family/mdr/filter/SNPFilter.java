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

public class SNPFilter implements SNPFilterInterface {

	protected MapFile mapData;
	protected int[] WSNP;

	protected int[] wseq = null;
	protected int[] bgseq = null;
	protected boolean filterFlag = false;
	ArrayList<SNP> snpList;
	protected HashSet<Integer> selectedSNPSet;
	protected HashSet<Integer> excludedSNPSet;
	protected HashSet<Integer> bgSNPSet;

	protected int[][] in_range = null;
	protected int[][] ex_range = null;

	public SNPFilter(MapFile mapData) {
		this.mapData = mapData;
		snpList = mapData.getMarkerList();
		selectedSNPSet = NewIt.newHashSet();
		excludedSNPSet = NewIt.newHashSet();
		bgSNPSet = NewIt.newHashSet();
	}

	public void Select() {
		
		if (Parameter.bgsnpFlag) {
			selectBackgroundSNP();
			filterFlag = true;
		}

		if (Parameter.inchrFlag || Parameter.exchrFlag) {
			selectChromosome();
			filterFlag = true;
		}

		if (Parameter.snpwindowFlag) {
			selectSNPWindow();
			filterFlag = true;
		}

		if (Parameter.snpPairFlag || Parameter.snpFlag) {
			selectSNPRange();
			selectSNPs();
			filterFlag = true;
		}

		if (Parameter.geneFlag) {
			selectGene();
			filterFlag = true;
		}

		if (Parameter.regionFlag) {
			selectSNPRegion();
			filterFlag = true;
		}

		if (filterFlag && selectedSNPSet.size() == 0 && excludedSNPSet.size() == 0 && bgSNPSet.size() == 0) {
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

		if (Parameter.inchrFlag) {
			for (int i = 0; i < snpList.size(); i++) {
				SNP snp = snpList.get(i);
				String chr = snp.getChromosome();
				int idx = ArrayUtils.indexOf(Parameter.in_chr, chr);
				if (idx >= 0) {
					includeSNP(i);
				}
			}
		}

		if (Parameter.exchrFlag) {
			for (int i = 0; i < snpList.size(); i++) {
				SNP snp = snpList.get(i);
				String chr = snp.getChromosome();
				int idx = ArrayUtils.indexOf(Parameter.ex_chr, chr);
				if (idx >= 0) {
					excludeSNP(i);
				}
			}
		}
	}

	private void selectSNPWindow() {

		ArrayList<SNP> snpwindowList = NewIt.newArrayList();
		ArrayList<Integer> snpwindowIndex = NewIt.newArrayList();

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
			System.err.println("could not identify the window size.");
		}
		for (int i = 0; i < snpList.size(); i++) {
			SNP snp = snpList.get(i);
			for (int j = 0; j < snpwindowList.size(); j++) {
				int d = snpwindowList.get(j).compareTo(snp);
				if (d != Integer.MAX_VALUE && d != Integer.MIN_VALUE) {
					int idx = snpwindowIndex.get(j).intValue();
					if (d > Parameter.snp_window[idx][0] && d < Parameter.snp_window[idx][1])
						includeSNP(i);
				}
			}
		}
	}
	
	private void selectGene() {
		
		ArrayList<ArrayList<Integer>> xsnps = NewIt.newArrayList();
		for (int i = 0; i < Parameter.gene.length; i++) {
			ArrayList<Integer> s = NewIt.newArrayList();
			xsnps.add(s);
		}
		for (int i = 0; i < snpList.size(); i++) {
			SNP snp = snpList.get(i);
			String chr = snp.getChromosome();
			int pos = snp.getPosition();
			int idx = ArrayUtils.indexOf(Parameter.gene_chr, chr);
			if (idx >= 0) {
				if ( pos >= (Parameter.gene_begin[idx] - Parameter.genewindow)*1000  && pos <= (Parameter.gene_end[idx] + Parameter.genewindow) *1000  ) {
					xsnps.get(idx).add(i);
					includeSNP(i);
				}
			}
		}		

		if (Parameter.snp2genefilesFlag) {
			for (int i = 0; i < xsnps.size(); i++) {
				StringBuffer sb = new StringBuffer(Parameter.gene[i] + ".gene");
				PrintStream PW = null;
				try {
					PW = new PrintStream(sb.toString());
				} catch (FileNotFoundException e) {
					e.printStackTrace();
				}

				ArrayList<Integer> s = xsnps.get(i);
				if (s.size() == 0) {
					continue;
				}
				System.err.println(s.size() + " snps selected with gene "
						+ Parameter.gene[i]);
				System.err.println("writing snps into " 
						+ Parameter.gene[i] + ".gene");
				Test.LOG.append(s.size() + " snps selected with gene "
						+ Parameter.gene[i] + "\n");
				Test.LOG.append("writing snps into " 
						+ Parameter.gene[i] + ".gene\n");
				for (int j = 0; j < s.size(); j++) {
					SNP snp = snpList.get(s.get(j));
					PW.println(snp.getName() + " " + snp.getChromosome() + " "
							+ snp.getPosition() +" " + Parameter.gene[i]);
				}
				PW.close();
			}
			Test.printLog();
			System.exit(1);
		} else {
			StringBuffer sb = new StringBuffer(Parameter.out + ".gene");
			PrintStream PW = null;
			try {
				PW = new PrintStream(sb.toString());
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
			for (int i = 0; i < xsnps.size(); i++) {
				ArrayList<Integer> s = xsnps.get(i);
				if (s.size() == 0) {
					continue;
				}
				System.err.println(s.size() + " snps selected with gene "
						+ Parameter.gene[i]);

				Test.LOG.append(s.size() + " snps selected with gene "
						+ Parameter.gene[i] + "\n");
				for (int j = 0; j < s.size(); j++) {
					SNP snp = snpList.get(s.get(j));
					PW.println(snp.getName() + " " + snp.getChromosome() + " "
							+ snp.getPosition() +" " + Parameter.gene[i]);
				}
			}
			System.err.println("writing snps into " 
					+ Parameter.out + ".gene");
			Test.LOG.append("writing snps into " 
					+ Parameter.out + ".gene\n");
			PW.close();
			if (Parameter.snp2genefileFlag) {
				Test.printLog();
				System.exit(1);
			}
		}

	}

	private void selectSNPRange() {

		if (Parameter.insnpPair != null) {
			in_range = new int[Parameter.insnpPair.length / 2][2];
			for (int i = 0; i < snpList.size(); i++) {
				SNP snp = snpList.get(i);
				String rs = snp.getName();
				int idx = ArrayUtils.indexOf(Parameter.insnpPair, rs);
				if (idx >= 0) {
					int dim1 = idx >> 1;
					int dim2 = (idx - (idx >> 1 << 1));
					in_range[dim1][dim2] = i;
					includeSNP(i);
				}
			}
			for (int i = 0; i < in_range.length; i++) {
				if (in_range[i][0] > in_range[i][1]) {
					in_range[i][0] = in_range[i][0]^in_range[i][1];
					in_range[i][1] = in_range[i][0]^in_range[i][1];
					in_range[i][0] = in_range[i][0]^in_range[i][1];
				}
			}
		}

		if (Parameter.exsnpPair != null) {
			ex_range = new int[Parameter.exsnpPair.length / 2][2];
			for (int i = 0; i < snpList.size(); i++) {
				SNP snp = snpList.get(i);
				String rs = snp.getName();
				int idx = ArrayUtils.indexOf(Parameter.exsnpPair, rs);
				if (idx >= 0) {
					int dim1 = idx >> 1;
					int dim2 = (idx - (idx >> 1 << 1));
					ex_range[dim1][dim2] = i;
					excludeSNP(i);
				}
			}
			for (int i = 0; i < ex_range.length; i++) {
				if (ex_range[i][0] > ex_range[i][1]) {
					ex_range[i][0] = ex_range[i][0]^ex_range[i][1];
					ex_range[i][1] = ex_range[i][0]^ex_range[i][1];
					ex_range[i][0] = ex_range[i][0]^ex_range[i][1];
				}
			}
		}

		for (int i = 0; i < snpList.size(); i++) {
			boolean flag = false;
			if (ex_range != null) {
				for (int j = 0; j < ex_range.length; j++) {
					if (i >= ex_range[j][0] && i <= ex_range[j][1]) {
						excludeSNP(i);
						flag = true;
						break;
					}
				}
			}
			if (!flag) {
				if (in_range != null) {
					for (int j = 0; j < in_range.length; j++) {
						if (i >= in_range[j][0] && i <= in_range[j][1]) {
							includeSNP(i);
						}
					}
				}
			}
		}
	}
	
	private void selectSNPRegion() {
		for(int i = 0; i < snpList.size(); i++) {
			SNP snp = snpList.get(i);
			String chr = snp.getChromosome();
			int pos = snp.getPosition();
			int idx = ArrayUtils.indexOf(Parameter.chr_reg, chr);
			if (idx >= 0) {
				if (pos >= Parameter.begin[idx] * 1000 && pos <= Parameter.end[idx] * 1000 ) {
					includeSNP(i);
				}
			}
		}
		return;
	}

	private void selectSNPs() {

		if (Parameter.includesnp != null) {
			for (int i = 0; i < snpList.size(); i++) {
				SNP snp = snpList.get(i);
				String rs = snp.getName();
				int idx = ArrayUtils.indexOf(Parameter.includesnp, rs);
				if (idx >= 0) {
					includeSNP(i);
				}
			}
		}

		if (Parameter.excludesnp != null) {
			for (int i = 0; i < snpList.size(); i++) {
				SNP snp = snpList.get(i);
				String rs = snp.getName();
				int idx = ArrayUtils.indexOf(Parameter.excludesnp, rs);
				if (idx >= 0) {
					excludeSNP(i);
				}
			}
		}

		return;
	}

	private void makeWSNPList() {

		if (selectedSNPSet.size() > 0) {
			if (selectedSNPSet.size() > bgSNPSet.size()) {
				WSNP = new int[selectedSNPSet.size()];
				int c = 0;
				for (Iterator<Integer> e = selectedSNPSet.iterator(); e.hasNext();) {
					Integer V = e.next();
					WSNP[c++] = V.intValue();
				}
			} else if (Parameter.order > bgSNPSet.size()) {
				WSNP = new int[snpList.size()];
				for (int i = 0; i < snpList.size(); i++) {
					WSNP[i] = i;
				}
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

		if (bgSNPSet.size() > 0) {
			bgseq = new int[bgSNPSet.size()];
			wseq = new int[WSNP.length - bgseq.length];
			int c1 = 0, c2 = 0;
			for (int i = 0; i < WSNP.length; i++) {
				if (bgSNPSet.contains(WSNP[i])) {
					bgseq[c1++] = i;
				} else {
					wseq[c2++] = i;
				}
			}
		} else {
			wseq = new int[WSNP.length];
			for (int i = 0; i < WSNP.length; i++) {
				wseq[i] = i;
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
