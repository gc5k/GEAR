package family.mdr.filter;

import java.util.HashSet;

import util.NewIt;

public class SNPFilterII {
	private SNPFilterInterface snpFilter;
	private HashSet<Integer> snpSet;
	private int[] wseq;
	private int[] bgSNP;

	public SNPFilterII (SNPFilterInterface snpFilter) {
		this.snpFilter = snpFilter;
		this.bgSNP = snpFilter.getBgSeq();

		snpSet = NewIt.newHashSet();
		wseq = snpFilter.getWSeq();
	}

	public void Filter() {
		
	}

	public int[] getwseq() {
		return wseq;
	}
	
	public int[] getBgSNP() {
		return bgSNP;
	}
	
	
}
