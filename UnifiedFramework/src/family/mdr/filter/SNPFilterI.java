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
	
	protected int[][] snps = null;
	
	public SNPFilterI(MapFile mapData) {
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

		ArrayList<HashSet<Integer>> chrSNPSet = NewIt.newArrayList();
		int L = Parameter.chr.length;
		for (int i = 0; i < L; i++) {
			HashSet<Integer> chrSet = NewIt.newHashSet(); 
			chrSNPSet.add(chrSet);
		}

		for (int i = 0; i < snpList.size(); i++) {
			SNP snp = snpList.get(i);
			String chr = snp.getChromosome();
			int idx = ArrayUtils.indexOf(Parameter.chr, chr);
			if (idx >= 0) {
				includeSNP(i);
				HashSet<Integer> chrSet = chrSNPSet.get(idx);
				chrSet.add(new Integer(i));
			}
		}

		chrSNP = new int[L][];
		for (int i = 0; i < L; i++) {
			HashSet<Integer> chrSet = chrSNPSet.get(i);
			if(chrSet.size() == 0) throw new IllegalArgumentException("could not find chromosome " + Parameter.chr[i]);
			snpArrays.add(chrSet);
			chrSNP[i] = new int[chrSet.size()];
			int c = 0;
			for (Iterator<Integer> e = chrSet.iterator(); e.hasNext(); ) {
				chrSNP[i][c++] = e.next().intValue();
			}
		}
	}

	private void selectSNPWindow() {
		ArrayList<SNP> snpwindowList = NewIt.newArrayList();
		ArrayList<Integer> snpwindowIndex = NewIt.newArrayList();
		ArrayList<HashSet<Integer>> snpwindowSet= NewIt.newArrayList();

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
			throw new IllegalArgumentException("could not find the snp windows specified.");
		}
		for (int i = 0; i < snpList.size(); i++) {
			SNP snp = snpList.get(i);
			for (int j = 0; j < snpwindowList.size(); j++) {
				int d = snpwindowList.get(j).compareTo(snp);
				if (d != Integer.MAX_VALUE && d != Integer.MIN_VALUE) {
					int idx = snpwindowIndex.get(j).intValue();
					if (d > Parameter.snp_window[idx][0] && d < Parameter.snp_window[idx][1]) {
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
				throw new IllegalArgumentException("could not find snps near " + Parameter.snpwindow[snpwindowIndex.get(i).intValue()]);
			}
			snpArrays.add(swSet);
			snpWin[i] = new int[snpwindowSet.get(i).size()];
			c = 0;
			for (Iterator<Integer> e = swSet.iterator(); e.hasNext(); ) {
				snpWin[i][c++] = e.next().intValue();
			}
		}
	}

	private void selectSNPRange() {

		ArrayList<HashSet<Integer>> inRangeSet = NewIt.newArrayList();
		for (int i = 0; i < Parameter.insnpPair.length / 2; i++) {
			HashSet<Integer> rSet = NewIt.newHashSet();
			inRangeSet.add(rSet);
		}

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

		for (int i = 0; i < snpList.size(); i++) {
			if (in_range != null) {
				for (int j = 0; j < in_range.length; j++) {
					if (i >= in_range[j][0] && i <= in_range[j][1]) {
						includeSNP(i);
						HashSet<Integer> rSet = inRangeSet.get(j);
						rSet.add(new Integer(i));
						break;
					}
				}
			}
		}
		
		for (int i = 0; i < inRangeSet.size(); i++) {
			HashSet<Integer> rSet = inRangeSet.get(i);
			if(rSet.size() == 0) throw new IllegalArgumentException("could not find the snp range " + in_range[i][0] + " " + in_range[i][1]);
			snpArrays.add(rSet); 
		}
	}

	private void selectSNPs() {

		snps = new int[Parameter.includesnp.length][1];
		Arrays.fill(snps, -1);
		if (Parameter.includesnp != null) {
			for (int i = 0; i < snpList.size(); i++) {
				SNP snp = snpList.get(i);
				String rs = snp.getName();
				int idx = ArrayUtils.indexOf(Parameter.includesnp, rs);
				if (idx >= 0) {
					includeSNP(i);
					snps[idx][0] = i;
				}
			}
		}

		for (int i = 0; i < snps.length; i++) {
			if (snps[i][0] == -1) {
				throw new IllegalArgumentException("could not find snp " + Parameter.includesnp[i]);
			}
			HashSet<Integer> snpSet = NewIt.newHashSet();
			snpSet.add(new Integer(snps[i][0]));
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
			if(Parameter.x) {
				throw new IllegalArgumentException("no snps selected for detecting interaction.");
			}
		}

		if (bgSNPSet.size() > 0) {
			bgseq = new int[bgSNPSet.size()];
			int c = 0;
			for (Iterator<Integer> e = bgSNPSet.iterator(); e.hasNext(); ) {
				int v = e.next().intValue();
				int idx = ArrayUtils.indexOf(WSNP, v);
				bgseq[c++] = v;
			}
			for (Iterator<Integer> e = bgSNPSet.iterator(); e.hasNext(); ) {
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
			L += e.size();
		}
		wseq = new int[L][];
		int c = 0;
		for (HashSet<Integer> e : snpArrays) {
			if (e.size() > 0) {
				wseq[c++] = new int[e.size()];
				int cc = 0;
				for (Integer I : e) {
					wseq[c][cc++] = I.intValue();
				}
				Arrays.sort(wseq[c]);
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
