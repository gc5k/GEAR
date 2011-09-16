package family.mdr.arsenal;

import java.math.BigInteger;

import admixture.parameter.Parameter;

public class ModelGeneratorII extends ModelGenerator {

	private int[][] interactionSeq;

	public ModelGeneratorII(int[][] interactionSeq, int[] bgseq) {
		super(bgseq);
		this.interactionSeq = interactionSeq;
	}

	@Override
	public void revup(int R) {

		comb = new int[interactionSeq.length];
		r = interactionSeq.length;
		a = new int[r];

		BigInteger nFact = BigInteger.ONE;
		for (int i = 0; i < r; i++) {
			nFact = nFact.multiply(new BigInteger(Integer.toString(interactionSeq[i].length)));
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
		while (a[i] == interactionSeq[i].length - 1) {
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
			comb[i + offset] = interactionSeq[i][s[i]];
			if (i != s.length - 1) {
				sb.append(interactionSeq[i][s[i]]);
				sb.append(MDRConstant.seperator);
			} else {
				sb.append(interactionSeq[i][s[i]]);
			}
		}
		return sb.toString();
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		int[][] seq = { { 1, 2, 3 }, { 4, 5, 6, 12 }, { 7, 8, 9, 11 } };
		ModelGenerator MG = new ModelGeneratorII(seq, null);
		MG.revup(seq.length);
		for (; MG.hasNext();) {
			String m = MG.next();
			System.out.println(m);
		}
	}

}
