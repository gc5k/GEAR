package family.pedigree.file;

import java.util.Arrays;
import java.util.HashMap;

import admixture.parameter.Parameter;

import util.NewIt;

/**
 * 
 * @author Guo-Bo Chen chenguobo@gmail.com Thanks Jelai Wang.
 */
public class SNP {
	private String chr = "";
	private String name = "";
	private float distance = 0;
	private int position = 0;
	private String minor;
	private char[] snp;

	private HashMap<String, String> genotypeHash = NewIt.newHashMap();

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

	public void setAllele(short[] freq) {

		StringBuffer sb = new StringBuffer();
		sb.append(snp[0]);
		sb.append(snp[0]);
		genotypeHash.put(new Integer(0).toString(), sb.toString());

		StringBuffer sb1 = new StringBuffer();
		if(snp[0] <= snp[1]) {
			sb1.append(snp[0]);
			sb1.append(snp[1]);
		} else {
			sb1.append(snp[1]);
			sb1.append(snp[0]);
		}
		genotypeHash.put((new Integer(1)).toString(), sb1.toString());

		StringBuffer sb3 = new StringBuffer();
		sb3.append(snp[1]);
		sb3.append(snp[1]);
		genotypeHash.put(new Integer(3).toString(), sb3.toString());

		StringBuffer sb4 = new StringBuffer();
		double d0 = freq[0];
		double d1 = freq[1];
		if(freq[0] < freq[1]) {
			sb4.append(snp[0] + "(" + String.format("%.2f", d0/(d0+d1)) + "), " + snp[1] + "(" + String.format("%.2f", d1/(d0+d1)) + ")");
		} else {
			sb4.append(snp[1] + "(" + String.format("%.2f", d1/(d0+d1)) + "), " + snp[0] + "(" + String.format("%.2f", d0/(d0+d1)) + ")");
		}
		minor = sb4.toString();
	}

	public void setAllele(char[] a, short[] freq) {
		if (a.length >= 3) {
			System.err.println("more than 2 alleles for " + name);
			System.exit(0);
		} else if (a[0] == Parameter.missing_allele.charAt(0)|| a[1] == Parameter.missing_allele.charAt(0)){
			System.err.println("more than 2 alleles for " + name);
			System.exit(0);
		}
		snp = new char[2];
		System.arraycopy(a, 0, snp, 0, 2);
		Arrays.sort(a);
		StringBuffer sb1 = new StringBuffer();
		sb1.append(a[0]);
		sb1.append(a[0]);
		genotypeHash.put((new Integer(0)).toString(), sb1.toString());
		
		StringBuffer sb2 = new StringBuffer();
		if (a[0]>=a[1]) {
			sb2.append(a[0]);
			sb2.append(a[1]);
		} else {
			sb2.append(a[1]);
			sb2.append(a[0]);
		}
		genotypeHash.put((new Integer(1)).toString(), sb2.toString());

		StringBuffer sb3 = new StringBuffer();
		sb3.append(a[1]);
		sb3.append(a[1]);
		genotypeHash.put((new Integer(3)).toString(), sb3.toString());
		
		StringBuffer sb = new StringBuffer();
		double d0 = freq[0];
		double d1 = freq[1];
		if(freq[0] < freq[1]) {
			sb.append(a[0] + "(" + String.format("%.2f", d0/(d0+d1)) + "), " + a[1] + "(" + String.format("%.2f", d1/(d0+d1)) + ")");
		} else {
			sb.append(a[1] + "(" + String.format("%.2f", d1/(d0+d1)) + "), " + a[0] + "(" + String.format("%.2f", d0/(d0+d1)) + ")");
		}
		minor = sb.toString();
	}

	public String getPolymorphism(String g) {
		return genotypeHash.get(g);
	}
	
	public char[] getSNP() {
		return snp;
	}

	public String toString() {
		StringBuffer sb = new StringBuffer();
		sb.append(name + " ");
		sb.append(chr + " ");
		sb.append(position + " ");
		if(minor != null) {
			sb.append(minor);
		}
		return sb.toString();
	}
}
