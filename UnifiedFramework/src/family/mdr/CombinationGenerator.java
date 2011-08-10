package family.mdr;

import java.util.Iterator;
import java.math.BigInteger;

import publicAccess.PublicData;

public class CombinationGenerator implements Iterator<String> {

	private int[] seq;
	private int[] a;
	private int[] comb;
	private String header = null;
	private int offset = 0;
	private int[] in_snp = null;
	private int[] ex_snp = null;
	private int n;
	private int r;
	private int len;
	private BigInteger numLeft;
	private BigInteger total;

	public CombinationGenerator(int n, int[] include_snp, int[] exclude_snp) {
		this.n = n;

		if (include_snp != null) {
			in_snp = include_snp;
			offset = in_snp.length;
		}

		if (exclude_snp != null) {
			ex_snp = exclude_snp;
		}

	}

	public void revup(int R) {

		r = R;
		len = n;
		comb = new int[R];
		if (in_snp != null) {
			len -= in_snp.length;
			r -= in_snp.length;
			System.arraycopy(in_snp, 0, comb, 0, offset);
			StringBuilder sb = new StringBuilder();
			for(int i = 0; i < in_snp.length; i++) {
				if(i != in_snp.length -1) {
					sb.append(in_snp[i] + PublicData.seperator);
				} else {
					sb.append(in_snp[i]);
				}
			}
			header = sb.toString();
		}
		
		if (ex_snp != null) {
			len -= ex_snp.length;
		}

		if (r < 0 || len < 0 || r > len) {
			System.err.println("impossible to draw " + r + " from " + len + " factors");
			throw new IllegalArgumentException();
		}

		a = new int[r];

		int idx = 0;
		seq = new int[len];
		for (int i = 0; i < n; i++) {
			boolean flag = true;
			if (in_snp != null) {
				for (int j = 0; j < in_snp.length; j++) {
					if (i == in_snp[j]) {
						flag = false;
						break;
					}
				}
			}
			if (!flag)
				continue;
			if (ex_snp != null) {
				for (int j = 0; j < ex_snp.length; j++) {
					if (i == ex_snp[j]) {
						flag = false;
						break;
					}
				}
			}
			if (flag) {
				seq[idx++] = i;
			}
		}

		BigInteger nFact = getFactorial(len);
		BigInteger rFact = getFactorial(r);
		BigInteger nminusrFact = getFactorial(len - r);
		total = nFact.divide(rFact.multiply(nminusrFact));
		reset();
	}

	public void reset() {
		for (int i = 0; i < a.length; i++) {
			a[i] = i;
		}
		numLeft = new BigInteger(total.toString());
	}

	private static BigInteger getFactorial(int n) {
		BigInteger fact = BigInteger.ONE;
		for (int i = n; i > 1; i--) {
			fact = fact.multiply(new BigInteger(Integer.toString(i)));
		}
		return fact;
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
		// return getString(a);
	}

	@Override
	public void remove() {
		throw new UnsupportedOperationException();
	}

	private String getComb(int[] s) {
		StringBuffer sb = new StringBuffer();
		if (header != null) {
			sb.append(header + PublicData.seperator);
		}
		for (int i = 0; i < s.length; i++) {
			comb[i + offset] = seq[s[i]];
			if(i != s.length - 1) {
				sb.append(seq[s[i]] + PublicData.seperator);
			} else {
				sb.append(seq[s[i]]);
			}
		}
		return sb.toString();

	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		int[] in = {0};
		int[] ex = {2};
		long t = System.currentTimeMillis();
		int len = 10;
		CombinationGenerator cg = new CombinationGenerator(len, in, ex);
		cg.revup(11);
		for (; cg.hasNext();) {
			System.out.println(cg.next());
		}
		System.out.println(System.currentTimeMillis() - t);
		System.out.println("done");
	}
}
