package family.pedigree.genotype;

import family.mdr.arsenal.MDRConstant;

/**
 * stores the genotypes of each individual. this class is not thread safe
 * (untested)
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class BPerson {

	protected String familyID;
	protected String personID;
	protected String momID;
	protected String dadID;
	protected int gender;
	protected String affectedStatus;
	protected int numMarkers;
	protected int genoLen;
	protected int[] alleles;
	protected final int intL = 16;
	protected final int shift = 4;
	protected final int mask = 15;

	public BPerson(int numMarkers) {
		this.numMarkers = numMarkers;
		if (numMarkers % intL == 0) {
			genoLen = numMarkers / intL;
		} else {
			genoLen = numMarkers / intL + 1;
		}
		alleles = new int[genoLen];
	}

	public BPerson(BPerson p) {
		familyID = p.getFamilyID();
		personID = p.getPersonID() + "ajhg2008";
		momID = p.getMomID();
		dadID = p.getDadID();
		affectedStatus = p.getAffectedStatus();
		this.numMarkers = p.getNumMarkers();
		if (numMarkers % intL == 0) {
			genoLen = numMarkers / intL;
		} else {
			genoLen = numMarkers / intL + 1;
		}
		alleles = new int[genoLen];
	}

	/**
	 * gets the family ID
	 * 
	 * @return The familyID for this person
	 */
	public String getFamilyID() {
		return familyID;
	}

	/**
	 * sets the family ID
	 * 
	 * @param familyID
	 */
	public void setFamilyID(String familyID) {
		this.familyID = familyID;
	}

	/**
	 * gets the Person ID
	 * 
	 * @return The personID for this person
	 */
	public String getPersonID() {
		return personID;
	}

	/**
	 * sets the person ID
	 * 
	 * @param personID
	 */
	public void setPersonID(String personID) {
		this.personID = personID;
	}

	/**
	 * gets the momID for this person
	 * 
	 * @return momID
	 */
	public String getMomID() {
		return momID;
	}

	/**
	 * sets the momid
	 * 
	 * @param momID
	 */
	public void setMomID(String momID) {
		this.momID = momID;
	}

	/**
	 * gets the dad ID for this person
	 * 
	 * @return dadID
	 */
	public String getDadID() {
		return dadID;
	}

	/**
	 * sets the dadID
	 * 
	 * @param dadID
	 */
	public void setDadID(String dadID) {
		this.dadID = dadID;
	}

	/**
	 * gets the gender for this person
	 * 
	 * @return gender
	 */
	public int getGender() {
		return gender;
	}

	/**
	 * sets the gender
	 * 
	 * @param gender
	 */
	public void setGender(int gender) {
		this.gender = gender;
	}

	/**
	 * gets the affected status for this person
	 * 
	 * @return affectedStatus
	 */
	public String getAffectedStatus() {
		return affectedStatus;
	}

	/**
	 * sets the affected status
	 * 
	 * @param affectedStatus
	 */
	public void setAffectedStatus(String affectedStatus) {
		this.affectedStatus = affectedStatus;
	}

	/**
	 * returns the number of markers for this person
	 * 
	 * @return integer count of markers
	 */
	public int getNumMarkers() {
		return numMarkers;
	}

	public void addMarker(boolean flag, int a1, int a2, int i) {
		int posByte = i >> shift;
		int posBite = (i & 0xf) << 1;

		if (flag) {
			if (a2 == a1) {// add 00 or 11
				int c = ((a1 << 1) + a2) << posBite;
				alleles[posByte] |= c;
			} else {// add 10
				int c = 2 << posBite;
				alleles[posByte] |= c;
			}
		} else {// add 01
			alleles[posByte] |= (1 << posBite);
		}
	}

	public void addByteGenotype(int g, int i) {
		int posByte = i >> shift;
		int posBite = (i & 0xf) << 1;

		alleles[posByte] |= g << posBite;
	}

	public void addAllMarker(byte[] g) {

		for (int i = 0; i < g.length; i++) {
			int posByte = i >> 2; // one int holds 4 bytes.
			int posBiteShift = (i - (i >> 2 << 2)) << 3;
			alleles[posByte] |= g[i] << posBiteShift;
		}

	}

	public String getGenotypeScoreString(int i) {
		int posByte = i >> shift;
		int posBite = (i & 0xf) << 1;
		int g = (alleles[posByte] >> (posBite)) & 3;
		if (g == 1) {// 01
			return MDRConstant.missingGenotype;
		} else {
			if (g == 2) {
				return Integer.toString(1);
			} else {
				return Integer.toString(g);
			}
		}
	}

	public int getAlleleArrayLength() {
		return genoLen;
	}

	public int[] getAlleleArray() {
		return alleles;
	}

	public String getBiAlleleGenotypeString(int i) {
		int posByte = i >> shift;
		int posBite = (i & 0xf) << 1;
		int g = (alleles[posByte] >> posBite) & 3;
		if (g == 1) {// 01
			return MDRConstant.missingGenotype;
		} else {
			StringBuffer sb = new StringBuffer();
			sb.append((alleles[posByte] >> (posBite + 1)) & 1);
			sb.append(alleles[posByte] >> posBite & 1);
			return sb.toString();
		}
	}

	public void setNonTransmittedGenotype(int index, String geno) {
		int a = Integer.parseInt(geno.substring(0, 1));
		int b = Integer.parseInt(geno.substring(1, 2));
		boolean flag = geno.compareTo(MDRConstant.missingGenotype) == 0 ? false : true;
		addMarker(flag, a, b, index);
	}
}
