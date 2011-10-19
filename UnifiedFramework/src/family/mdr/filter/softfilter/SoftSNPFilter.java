package family.mdr.filter.softfilter;

import java.util.HashSet;

import admixture.parameter.Parameter;

import family.mdr.filter.SNPFilterInterface;
import family.popstat.AlleleFrequency;

import util.NewIt;

public class SoftSNPFilter {
	
	private HashSet<Integer> qualifiedSNPSet;

	private AlleleFrequency al;
	private int[] wseq;
	private int[][] wseq2;
	private int[] bgSNP;

	public SoftSNPFilter (SNPFilterInterface snpFilter, AlleleFrequency al) {

		this.al = al;
		this.bgSNP = snpFilter.getBgSeq();

		qualifiedSNPSet = NewIt.newHashSet();
		wseq = snpFilter.getWSeq();
		wseq2 = snpFilter.getWSeq2();
		
	}

	public void Filter() {
		if (wseq != null) {
			double [][] allelefreq = al.getAlleleFrequency();
			for (int i = 0; i < wseq.length; i++) {
				int idx = wseq[i];
				boolean flag = true;
				if (Parameter.genoFlag) {
					if (allelefreq[idx][2] <= Parameter.geno) {
						qualifiedSNPSet.add(new Integer(idx));
					} else {
						flag = false;
					}
				}
				if (Parameter.mafFlag) {
					double f = allelefreq[idx][0] < allelefreq[idx][1] ? allelefreq[idx][0]:allelefreq[idx][1];
					if (f > Parameter.maf && f != 0) {
						qualifiedSNPSet.add(new Integer(idx));
					} else {
						flag = false;
					}
				}
				if (flag) {
					qualifiedSNPSet.add(new Integer(idx));
				}
			}
			
			wseq = new int[qualifiedSNPSet.size()];
			int c = 0;
			for (Integer I : qualifiedSNPSet) {
				wseq[c++] = I.intValue();
			}
		}
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
