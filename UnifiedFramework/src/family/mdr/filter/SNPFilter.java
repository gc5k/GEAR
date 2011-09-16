package family.mdr.filter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;

import org.apache.commons.lang3.ArrayUtils;

import admixture.parameter.Parameter;
import family.pedigree.file.MapFile;
import family.pedigree.file.SNP;

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
		}

		if (Parameter.chrFlag) {
			selectChromosome();
		}

		if (Parameter.snpwindowFlag) {
			selectSNPWindow();
		}

		if (Parameter.snpPairFlag || Parameter.snpFlag) {
			selectSNPRange();
			selectSNPs();
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

		for (int i = 0; i < snpList.size(); i++) {
			SNP snp = snpList.get(i);
			String chr = snp.getChromosome();
			int idx = ArrayUtils.indexOf(Parameter.chr, chr);
			if (idx >= 0) {
				includeSNP(i);
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
			throw new IllegalArgumentException("could not identify the window size.");
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

		if (selectedSNPSet.size() > 0 ) {
			WSNP = new int[selectedSNPSet.size()];
			int c = 0;
			for (Iterator<Integer> e = selectedSNPSet.iterator(); e.hasNext();) {
				Integer V = e.next();
				WSNP[c++] = V.intValue();
			}
			Arrays.sort(WSNP);
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
				if (bgSNPSet.contains(new Integer(WSNP[i]))) {
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
		Integer I = new Integer(i);
		if (selectedSNPSet.contains(I)) {
			return true;
		} else {
			selectedSNPSet.add(I);
			return false;
		}
	}

	private boolean excludeSNP(int i) {
		Integer I = new Integer(i);
		if (excludedSNPSet.contains(I)) {
			return true;
		} else {
			excludedSNPSet.add(I);
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
		// TODO Auto-generated method stub
		return null;
	}
}
