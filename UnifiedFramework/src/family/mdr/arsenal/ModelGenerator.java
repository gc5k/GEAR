package family.mdr.arsenal;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.math.BigInteger;

import admixture.parameter.Parameter;
import test.Test;
import util.NewIt;

public class ModelGenerator implements Iterator<String> {

	protected int[] seq;
	protected int[] a;
	protected int[] comb;
	protected String header = null;
	protected int offset = 0;

	protected int[] bgSNP;
	protected int r;
	protected int len;
	protected BigInteger numLeft;
	protected BigInteger total;

	public ModelGenerator(int[] bgSNP) {
		this.bgSNP = bgSNP;
	}

	public ModelGenerator(int[] seq, int[] bgSNP) {
		this.seq = seq;
		this.bgSNP = bgSNP;
	}

	public void revup(int R) {

		if (bgSNP != null) {
			if (R > bgSNP.length + seq.length) {

				System.err.println("the number of qualified snps were fewer than the order of interaction going to detect.");
				Test.LOG.append("the number of qualified snps were fewer than the order of interaction going to detect.\n");
				Test.printLog();
				System.exit(0);
			}
			if (R < bgSNP.length) {
				R = bgSNP.length;
			}

		} else {
			if (R > seq.length) {
				System.err.println("the number of qualified snps were fewer than the order of interaction going to detect.");
				Test.LOG.append("the number of qualified snps were fewer than the order of interaction going to detect.\n");
				Test.printLog();
				System.exit(0);				
			}
		}
		len = seq.length;
		comb = new int[R];

		if (bgSNP != null) {
			r = R - bgSNP.length;
			System.arraycopy(bgSNP, 0, comb, 0, bgSNP.length);
			StringBuilder sb = new StringBuilder();
			for (int i = 0; i < bgSNP.length; i++) {
				if (i != bgSNP.length - 1) {
					sb.append(bgSNP[i]);
					sb.append(MDRConstant.seperator);
				} else {
					sb.append(bgSNP[i]);
				}
			}
			header = sb.toString();
		} else {
			r = R;
		}

		a = new int[r];

		BigInteger nFact = BigInteger.ONE;
		for (int i = 0; i < r; i++) {
			nFact = nFact.multiply(new BigInteger(Integer.toString(len - i)));
		}
		for (int i = 1; i <= r; i++) {
			nFact = nFact.divide(new BigInteger(Integer.toString(i)));
		}

		total = nFact;
		reset();
	}

	public void reset() {
		for (int i = 0; i < a.length; i++) {
			a[i] = i;
		}
		numLeft = new BigInteger(total.toString());
		if (Parameter.sliceN > 1) {
			BigInteger size = total.divide(new BigInteger(Integer.toString(Parameter.sliceN)));
			if (size.compareTo(BigInteger.ZERO) == 0) {
				System.err.println("Given " + total + " interactions, " + "--slice option made " + Parameter.sliceN + " slices, which " +
						"were too many.\nGMDR quit.");
				Test.LOG.append("Given " + total + " interactions, " + "--slice option made " + Parameter.sliceN + " slices, which " +
						"were too many.\nGMDR quit\n");
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

	protected void shift(BigInteger stop) {
		while (stop.compareTo(BigInteger.ZERO) != 0) {
			inner_next();
			stop = stop.subtract(BigInteger.ONE);
		}
	}

	protected void inner_next() {
		int i = r - 1;
		while (a[i] == len - r + i) {
			i--;
		}
		a[i] = a[i] + 1;
		for (int j = i + 1; j < r; j++) {
			a[j] = a[i] + j - i;
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
		while (a[i] == len - r + i) {
			i--;
		}
		a[i] = a[i] + 1;
		for (int j = i + 1; j < r; j++) {
			a[j] = a[i] + j - i;
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
			comb[i + offset] = seq[s[i]];
			if (i != s.length - 1) {
				sb.append(seq[s[i]]);
				sb.append(MDRConstant.seperator);
			} else {
				sb.append(seq[s[i]]);
			}
		}
		return sb.toString();
	}

	public BigInteger getTotal() {
		return total;
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		System.out.println("CombinationGeneratorII");
		int[] in = { 2 };
		HashSet<Integer> S = NewIt.newHashSet();
		for (int i = 0; i < 100; i++) {
			S.add(new Integer(i));
		}
		int[] seq = new int[S.size()];

		int c = 0;
		for (Iterator<Integer> e = S.iterator(); e.hasNext();) {
			Integer I = e.next();
			seq[c++] = I.intValue();
		}
		Arrays.sort(seq);

		ModelGenerator cg = new ModelGenerator(seq, in);
		System.out.println(System.currentTimeMillis());
		long t1 = System.currentTimeMillis();
		cg.revup(2);
		long t2 = System.currentTimeMillis();
		System.out.println("start:" + (t2 - t1));
		for (; cg.hasNext();) {
			String g = cg.next();
			System.out.println(g);
		}
		System.out.println(System.currentTimeMillis() - t2);
		System.out.println("done");
	}
}
