package family.mdr.filter.softfilter;

import java.util.HashSet;

import family.mdr.filter.SNPFilterInterface;

import util.NewIt;

public class SoftSNPFilter {
	private SNPFilterInterface snpFilter;
	private HashSet<Integer> snpSet;
	private int[] wseq;
	private int[][] wseq2;
	private int[] bgSNP;

	public SoftSNPFilter (SNPFilterInterface snpFilter) {
		this.snpFilter = snpFilter;
		this.bgSNP = snpFilter.getBgSeq();

		snpSet = NewIt.newHashSet();
		wseq = snpFilter.getWSeq();
		wseq2 = snpFilter.getWSeq2();
	}

	public void Filter() {
		
	}

	public int[] getWSeq() {
		return wseq;
	}
	
	public int[][] getWSeq2() {
		return wseq2;
	}
	public int[] getBgSeq() {
		return bgSNP;
	}
	
	
}
