package family.pedigree.file;

import test.Test;
import admixture.parameter.Parameter;

/**
 * 
 * @author Guo-Bo Chen chenguobo@gmail.com
 * Thanks Jelai Wang.
 */
public class SNP implements Comparable<SNP>{
	private String chr = "";
	private String name = "";
	private float distance = 0;
	private int position = 0;
	private String minor;
	private char[] snp;
	private double[] freq;

	public SNP(String n) {
		name = n;
	}

	public SNP(String c, String n, float d, int p) {
		chr = c;
		name = n;
		distance = d;
		position = p;
	}

	public SNP(String c, String n, float d, int p, char a1, char a2) {
		chr = c;
		name = n;
		distance = d;
		position = p;
		snp = new char[2];
		snp[0] = a1;
		snp[1] = a2;
	}

	/**
	 * Returns the name of the SNP.
	 */
	public String getName() {
		return name;
	}

	/**
	 * Returns the chromosome where the SNP is located.
	 */
	public String getChromosome() {
		return chr;
	}

	/**
	 * Returns the physical position on the chromosome where the SNP is located.
	 */
	public int getPosition() {
		return position;
	}

	/**
	 * Returns the genetic distance on the chromosome where the SNP is located.
	 */
	public float getDistance() {
		return distance;
	}

	public void setAllele(char[] a, short[] freq) {
		if (a.length >= 3) {
			System.err.println("more than 2 alleles for " + name);
			Test.LOG.append("more than 2 alleles for " + name + ".\n");
			Test.printLog();
			System.exit(0);
		} else if (a[0] == Parameter.missing_allele.charAt(0)
				|| a[1] == Parameter.missing_allele.charAt(0)) {
//			System.err.println("more than 2 alleles for " + name);
//			System.exit(0);
		}
		snp = new char[2];
		System.arraycopy(a, 0, snp, 0, 2);

		this.freq = new double[freq.length];
		System.arraycopy(freq, 0, this.freq, 0, freq.length);
	}

	public void setAllele(double[] freq) {
		this.freq = new double[freq.length];
		System.arraycopy(freq, 0, this.freq, 0, freq.length);
//		System.err.println(name + " " + chr + " " + position + " " + freq[0] + " " + freq[1] + " " + snp[0] + " " + snp[1]);
	}

	public void setAllelePolymorphism(char[] a) {
		snp = new char[2];
		System.arraycopy(a, 0, snp, 0, 2);
	}

	public String getPolymorphism(String g) {
		StringBuffer sb = new StringBuffer();
		switch (Integer.parseInt(g)) {
		case 0:
			sb.append(snp[0]);
			sb.append(snp[0]);
			break;
		case 1:
			sb.append(snp[0]);
			sb.append(snp[1]);
			break;
		case 2:
			sb.append(snp[1]);
			sb.append(snp[0]);
		case 3:
			sb.append(snp[1]);
			sb.append(snp[1]);
		}
		return sb.toString();
	}

	public char[] getSNP() {
		return snp;
	}

	public String toString() {
		StringBuffer sb = new StringBuffer();
		sb.append(name + ", ");
		sb.append(chr + ", ");
		sb.append(position + ", ");
		double d0 = 0;
		double d1 = 0;
		if(freq.length>1) {
			d0 = freq[0];
			d1 = freq[1];
		} else {
			d0 = freq[0];
		}
		if (d0 < d1) {
			sb.append(snp[0] + ", " + String.format("%.3f", d0 / (d0 + d1)) + ", ");
//					+ "), " + snp[1] + "("
//					+ String.format("%.2f", d1 / (d0 + d1)) + ")");
		} else {
			sb.append(snp[1] + ", " + String.format("%.3f", d1 / (d0 + d1)) + ", ");
//					+ "), " + snp[0] + "("
//					+ String.format("%.2f", d0 / (d0 + d1)) + ")");
		}
		if (minor != null) {
			sb.append(minor);
		}
		return sb.toString();
	}

	@Override
	public int compareTo(SNP snp) {
		if(snp.chr.compareTo(this.chr) == 0) {
			return snp.getPosition() - position;
		} else if (snp.chr.compareTo(this.chr) < 0){
			return Integer.MIN_VALUE;
		} else {
			return Integer.MAX_VALUE;
		}
	}
}
