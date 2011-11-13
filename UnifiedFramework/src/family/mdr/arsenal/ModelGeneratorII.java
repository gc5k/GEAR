package family.mdr.arsenal;

import java.math.BigInteger;

import test.Test;

import admixture.parameter.Parameter;

public class ModelGeneratorII extends ModelGenerator {

	private int[][] wseq2;

	public ModelGeneratorII(int[][] wseq2, int[] bgseq) {
		super(bgseq);
		this.wseq2 = wseq2;
		if (bgseq != null) {
			Parameter.order = wseq2.length + bgseq.length;
		} else {
			Parameter.order = wseq2.length;
		}
	}

	@Override
	public void revup(int R) {
		int L = 0;
		if(bgSNP != null) {
			L += bgSNP.length;
		}
		comb = new int[wseq2.length + L];
		r = wseq2.length;
		a = new int[r];
		if (bgSNP != null) {
			r = R - bgSNP.length;
			System.arraycopy(bgSNP, 0, comb, 0, bgSNP.length);
			StringBuilder sb = new StringBuilder();
			for(int i = 0; i < bgSNP.length; i++) {
				if(i != bgSNP.length -1) {
					sb.append(bgSNP[i]);
					sb.append(MDRConstant.seperator);
				} else {
					sb.append(bgSNP[i]);
				}
			}
			header = sb.toString();
		}
		BigInteger nFact = BigInteger.ONE;
		for (int i = 0; i < r; i++) {
			nFact = nFact.multiply(new BigInteger(Integer.toString(wseq2[i].length)));
		}

		total = nFact;
		reset();
	}

	@Override
	public void reset() {
		for (int i = 0; i < a.length; i++) {
			a[i] = 0;
		}
		numLeft = new BigInteger(total.toString());
		if (Parameter.sliceN > 1) {
			BigInteger size = total.divide(new BigInteger(Integer.toString(Parameter.sliceN)));
			if (size.compareTo(BigInteger.ZERO) == 0) {
				System.err.println("Given " + total + " interactions, " + "--slice option made " + Parameter.sliceN + " slices, which " +
						"were too many.\nGMDR quit.");
				Test.LOG.append("Given " + total + " interactions, " + "--slice option made " + Parameter.sliceN + " slices, which " +
						"were too many.\nGMDR quit.\n");
				Test.printLog();
				System.exit(1);
			}
			BigInteger stop = size.multiply(new BigInteger(Integer.toString(Parameter.slice - 1)));
			shift(stop);
			if (Parameter.slice != Parameter.sliceN) {
				numLeft = size;
			} else {
				BigInteger s = size.multiply(new BigInteger(Integer.toString(Parameter.sliceN - 1)));
				numLeft = total.subtract(s);
				numLeft = numLeft.subtract(BigInteger.ONE);
			}
		}
	}

	@Override
	public boolean hasNext() {
		return numLeft.compareTo(BigInteger.ZERO) == 1;
	}

	@Override
	public String next() {

		if (numLeft.equals(total)) {
			numLeft = numLeft.subtract(BigInteger.ONE);
			return getComb(a);
		}

		int i = r - 1;
		while (a[i] == wseq2[i].length - 1) {
			i--;
		}
		a[i] = a[i] + 1;
		for (int j = i + 1; j < r; j++) {
			a[j] = 0;
		}

		numLeft = numLeft.subtract(BigInteger.ONE);

		return getComb(a);
	}

	@Override
	public void remove() {
		throw new UnsupportedOperationException();
	}

	private String getComb(int[] s) {
		StringBuffer sb = new StringBuffer();
		if (header != null) {
			sb.append(header);
			sb.append(MDRConstant.seperator);
		}
		for (int i = 0; i < s.length; i++) {
			comb[i + offset] = wseq2[i][s[i]];
			if (i != s.length - 1) {
				sb.append(wseq2[i][s[i]]);
				sb.append(MDRConstant.seperator);
			} else {
				sb.append(wseq2[i][s[i]]);
			}
		}
		return sb.toString();
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		int[][] seq = { { 1, 2, 3 }, { 4, 5, 6, 12 }, { 7, 8, 9, 11 } };
		ModelGenerator MG = new ModelGeneratorII(seq, null);
		MG.revup(seq.length);
		for (; MG.hasNext();) {
			String m = MG.next();
			System.out.println(m);
		}
	}

}
