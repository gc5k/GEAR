package family.pedigree.genotype;

import java.util.ArrayList;

import util.NewIt;

/**
 * stores the genotypes of each individual. this class is not thread safe (untested)
 * 
 * @author Guo-Bo Chen
 */
public class Person {

    private String familyID;
    private String personID;
    private String momID;
    private String dadID;
    private int gender;
    private int affectedStatus;
//    private ArrayList<String> genotype;
    private String reasonImAxed;
    // private Vector markers;
    // private byte[] alleles1;
    // private byte[] alleles2;
    private byte[][] alleles;
    private double numGoodMarkers;
    private boolean[] zeroed;
    // this is used to keep track of the index of the last marker added
    private int currMarker;
    public final static int FEMALE = 2;
    public final static int MALE = 1;
    public final static int AFFACTED = 2;
    public final static int UNAFFACTED = 1;
    public final static String DATA_MISSING = "0";

    public Person(int numMarkers) {
        alleles = new byte[2][numMarkers];
        this.zeroed = new boolean[numMarkers];
        this.currMarker = 0;
//        this.genotype = NewIt.newArrayList();
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

    public void setMarkers(byte[] ma, byte[] mb) {
        alleles[0] = ma;
        alleles[1] = mb;
    }

    /**
     * returns the number of markers for this person
     * 
     * @return integer count of markers
     */
    public int getNumMarkers() {
        return this.alleles[0].length;
    }

    public byte getAllele(int location, int index) {
        return alleles[index][location];
    }

    public void addMarker(byte markera, byte markerb) {
        alleles[0][currMarker] = markera;
        alleles[1][currMarker] = markerb;
        zeroed[currMarker] = false;
        currMarker++;
        if (!(markera == 0 || markerb == 0)) {
            numGoodMarkers++;
        }
    }

    /**
     * checks to see if a marker has been zeroed out
     * 
     * @param location
     *            - which marker to check
     * @return true if marker is zeroed, false otherwise
     */
    public boolean getZeroed(int location) {
        // return ((Boolean)zeroed.get(location)).booleanValue();
        return zeroed[location];
    }

    /**
     * sets the bit that this marker has been zeroed out for this indiv (e.g. because it has a mendel error)
     * 
     * @param i
     *            - marker to be zeroed
     */
    public void zeroOutMarker(int i) {
        // this.zeroed.set(i, new Boolean(true));
        this.zeroed[i] = true;
    }

    public String getReasonImAxed() {
        return reasonImAxed;
    }

    public void setReasonImAxed(String reasonImAxed) {
        this.reasonImAxed = reasonImAxed;
    }

    public double getGenoPC() {
        return numGoodMarkers / alleles[0].length;
    }

    public boolean[] getZeroedArray() {
        return zeroed;
    }

    public void setZeroedArray(boolean[] z) {
        zeroed = z;
    }

//    public void setGenotype(ArrayList<String> marker) {
//        genotype = marker;
//    }
//
    public void setGenotype(int index, String geno) {
    	alleles[0][index] = Byte.parseByte(geno.substring(0, 1));
    	alleles[1][index] = Byte.parseByte(geno.substring(1, 2));
    }
//
//    public ArrayList<String> getGenotype() {
//        return genotype;
//    }
//
    public ArrayList<String> getGenotype(int[] subsetMarker) {
    		ArrayList<String> sub = NewIt.newArrayList();
    		for (int i = 0; i < subsetMarker.length; i++) {
    			StringBuffer sb = new StringBuffer();
    			sb.append(alleles[0][subsetMarker[i]]);
    			sb.append(alleles[1][subsetMarker[i]]);
    			sub.add(sb.toString());
    		}
    		return sub; 	
    }

    public String getGenotype(int index) {
    	StringBuilder sb = new StringBuilder();
    	sb.append(alleles[0][index]);
    	sb.append(alleles[1][index]);
        return sb.toString();
    }
    
    public byte[][] getGenotype() {
    	return alleles;
    }
}
