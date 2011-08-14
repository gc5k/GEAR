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
	private String chr;
	private String name;
	private float distance;
	private int position;

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

	public void setAllele(String[] a) {
		if (a.length >= 3) {
			System.err.println("more than 2 alleles for " + name);
			System.exit(0);
		} else if (a[0].compareTo(Parameter.missing_allele)== 0|| a[1].compareTo(Parameter.missing_allele) == 0){
			System.err.println("more than 2 alleles for " + name);
			System.exit(0);
		}
		Arrays.sort(a);
		StringBuffer sb1 = new StringBuffer();
		sb1.append(a[0]);
		sb1.append(a[0]);
		genotypeHash.put((new Integer(0)).toString(), sb1.toString());
		
		StringBuffer sb2 = new StringBuffer();
		if (a[0].compareTo(a[1]) > 0) {
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
		genotypeHash.put((new Integer(2)).toString(), sb3.toString());

	}

	public String getPolymorphism(String g) {
		return genotypeHash.get(g);
	}

	public String toString() {
		StringBuffer sb = new StringBuffer();
		sb.append(name + " ");
		sb.append(chr + " ");
		sb.append(position);
		return sb.toString();
	}
}
