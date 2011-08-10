package family.pedigree.file;

/**
 * 
 * @author Guo-Bo Chen chenguobo@gmail.com
 * Thanks Jelai Wang.
 */
public class SNP {
	private String chr;
	private String name;
	private float distance;
	private int position;
	
	public SNP (String n) {
		name = n;
	}
	
	public SNP (String c, String n, float d, int p) {
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
	
	public String toString() {
		StringBuffer sb = new StringBuffer();
		sb.append(name + " ");
		sb.append(chr + " ");
		sb.append(position);
		return sb.toString();
	}
}
