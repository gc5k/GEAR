package family.pedigree.file;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Set;

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
	private String[] allele = null;

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

	public void setAllele(Set<String> a) {
		if (a.size() >= 3) {
			System.err.println("more than 2 alleles for " + name);
			System.exit(0);
		} else if (a.size() == 1){
			System.err.println("more than 2 alleles for " + name);
			System.exit(0);
		} else {
			allele = (String[]) a.toArray(new String[a.size()]);
		}
		Arrays.sort(allele);
		byte[] genotype = new byte[2];
		for (int i = 0; i < 2; i++) {
			if (allele[i].equalsIgnoreCase("A")) {
				genotype[i] = 1;
			} else if (allele[i].equalsIgnoreCase("C")) {
				genotype[i] = 2;
			} else if (allele[i].equalsIgnoreCase("G")) {
				genotype[i] = 3;
			} else if (allele[i].equalsIgnoreCase("T")) {
				genotype[i] = 4;
			}
		}
		StringBuffer sb1 = new StringBuffer();
		sb1.append(allele[0]);
		sb1.append(allele[0]);
		genotypeHash.put((new Integer(genotype[0] + genotype[0])).toString(), sb1.toString());
		
		StringBuffer sb2 = new StringBuffer();
		sb2.append(allele[0]);
		sb2.append(allele[1]);
		genotypeHash.put((new Integer(genotype[0] + genotype[1])).toString(), sb2.toString());

		StringBuffer sb3 = new StringBuffer();
		sb3.append(allele[1]);
		sb3.append(allele[1]);
		genotypeHash.put((new Integer(genotype[1] + genotype[1])).toString(), sb3.toString());

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
