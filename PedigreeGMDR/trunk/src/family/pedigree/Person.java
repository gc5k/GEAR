package family.pedigree;

import java.util.ArrayList;

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
    private ArrayList genotype;
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
        this.genotype = new ArrayList();
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

    public void setGenotype(ArrayList marker) {
        genotype = marker;
    }

    public void setGenotype(int index, String geno) {
        genotype.set(index, geno);
    }

    public ArrayList getGenotype() {
        return genotype;
    }

    public String getGenotype(int index) {
        return (String) genotype.get(index);
    }
}
