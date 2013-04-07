package family.qc.colqc;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;

import org.apache.commons.lang3.ArrayUtils;

import family.pedigree.file.MapFile;
import family.pedigree.file.SNP;
import gear.Parameter;
import gear.util.Logger;
import gear.util.NewIt;

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
		makeWSNPList();
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
		} else  if (Parameter.INSTANCE.transFlag) {
			Logger.printUserError("No SNP is selected.");
			System.exit(1);
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
			Logger.printUserError("No SNP is selected.");
			System.exit(1);
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
