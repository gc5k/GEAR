package family.mdr;

import java.util.Arrays;
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
	private int[] outlier = null;
	private int n;
	private int r;
	private int len;
	private BigInteger numLeft;
	private BigInteger total;

	public CombinationGenerator(int n, int[] include_snp, int[] exclude_snp) {
		this.n = n;
		int L = 0;
		if (include_snp != null) {
			in_snp = include_snp;
			offset = in_snp.length;
			L += in_snp.length;
		}

		if (exclude_snp != null) {
			ex_snp = exclude_snp;
			L += ex_snp.length;
		}
		if (L > 0) {
			outlier = new int[L];
			if(in_snp != null && ex_snp != null) {
				System.arraycopy(in_snp, 0, outlier, 0, in_snp.length);
				System.arraycopy(ex_snp, 0, outlier, in_snp.length, ex_snp.length);
			} else if(in_snp != null && ex_snp == null) {
				System.arraycopy(in_snp, 0, outlier, 0, in_snp.length);
			} else if(in_snp == null && ex_snp != null) {
				System.arraycopy(ex_snp, 0, outlier, 0, ex_snp.length);
			}
			Arrays.sort(outlier);
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

		if (r < 0 || len < 0 || r > len) {
			System.err.println("impossible to draw " + r + " from " + len + " factors");
			throw new IllegalArgumentException();
		}
		
		a = new int[r];
		
		if (outlier != null) {
			seq = new int[len-outlier.length];
		} else {
			seq = new int[len];
		}
		int c = 0;
		if(outlier != null) {
			for (int i = 0; i < len; i++) {
				int idx = Arrays.binarySearch(outlier, i);
				if (idx >= 0 ) seq[c++] = i;
			}
		} else {
			for (int i = 0; i < len; i++) {
				seq[i] = i;
			}
		}

		BigInteger nFact = BigInteger.ONE;
		for (int i = 0; i < r; i++) {
			nFact = nFact.multiply(new BigInteger(Integer.toString(len-i)));
		}
		for (int i = 1; i <= r; i++) {
			nFact = nFact.divide(new BigInteger(Integer.toBinaryString(i)));
		}

		total = nFact;
		reset();
	}

	public void reset() {
		for (int i = 0; i < a.length; i++) {
			a[i] = i;
		}
		numLeft = new BigInteger(total.toString());
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
			sb.append(PublicData.seperator);
		}
		for (int i = 0; i < s.length; i++) {
			comb[i + offset] = seq[s[i]];
			if(i != s.length - 1) {
				sb.append(seq[s[i]]);
				sb.append(PublicData.seperator);
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
		System.out.println("CombinationGeneratorII");
		int[] in = null;
		int[] ex = null;
		long t = System.currentTimeMillis();
		int len = 1000000;
		CombinationGenerator cg = new CombinationGenerator(len, in, ex);
		System.out.println(System.currentTimeMillis());
		long t1 = System.currentTimeMillis();
		cg.revup(1);
		long t2 = System.currentTimeMillis();
		System.out.println("start:" + (t2 - t1));
		for (; cg.hasNext();) {
			String g = cg.next();
//			System.out.println(g);
		}
		System.out.println(System.currentTimeMillis() - t2);
		System.out.println("done");
	}
}
