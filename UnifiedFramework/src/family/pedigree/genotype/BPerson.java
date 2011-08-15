package family.pedigree.genotype;

import family.mdr.MDRConstant;

/**
 * stores the genotypes of each individual. this class is not thread safe
 * (untested)
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class BPerson {

	private String familyID;
	private String personID;
	private String momID;
	private String dadID;
	private int gender;
	private int affectedStatus;
	private int numMarkers;
	private int genoLen;
	private int[][] alleles;
	private int intL = 32;
	private int shift = 5;

	public BPerson(int numMarkers) {
		this.numMarkers = numMarkers;
		if (numMarkers % intL == 0) {
			genoLen = numMarkers / intL;
		} else {
			genoLen = numMarkers / intL + 1;
		}
		alleles = new int[3][genoLen];
		// this.genotype = NewIt.newArrayList();
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
	public int getAffectedStatus() {
		return affectedStatus;
	}

	/**
	 * sets the affected status
	 * 
	 * @param affectedStatus
	 */
	public void setAffectedStatus(int affectedStatus) {
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

	public int isGenotype(int idx) {
		int posByte = idx >> shift;
		int posBite = idx - (idx >> shift << shift);
		int f = alleles[0][posByte];
		return 1 & f >> posBite;
	}

	public void addMarker(boolean flag, int a, int b, int idx) {
		int posByte = idx >> shift;
		int posBite = idx - (idx >> shift << shift);
		int f = alleles[0][posByte];
		int m = flag ? 1 : 0;
		alleles[0][posByte] = f | (m << posBite);

		int a1 = alleles[1][posByte];
		alleles[1][posByte] = a1 | (a << posBite);

		int a2 = alleles[2][posByte];
		alleles[2][posByte] = a2 | (b << posBite);
	}

	public String getGenotypeString(int i) {
		StringBuffer sb = new StringBuffer();
		int posByte = i >> shift;
		int posBite = i - (i >> shift << shift);
		if (alleles[0][posByte] >> posBite == 0) {
			sb.append(MDRConstant.missingGenotype);
		} else {
			sb.append((1 & alleles[1][posByte] >> posBite) + (1 & alleles[2][posByte] >> posBite));
		}
		return sb.toString();
	}
}
