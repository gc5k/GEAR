package family.pedigree.genotype;

import java.util.ArrayList;

import util.NewIt;

/**
 * stores the genotypes of each individual, and, furthermore, store the nontransmitted genotype of this individual. this
 * class is not thread safe (untested)
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class PseudoPerson {

    private String familyID;
    private String pseudopersonID;
    private String momID;
    private String dadID;
    private int gender;
    private int affectedStatus;
    private byte[][] alleles;
    private int numMarkers;

    public PseudoPerson() {

    }

    public void MakeSpace(int n) {
    	numMarkers = n;
    	alleles = new byte[2][n];
    }

    public void addMarker(int i, String genotype) {
        alleles[0][i] = Byte.parseByte(genotype.substring(0,1));
        alleles[1][i] = Byte.parseByte(genotype.substring(1, 2));
    }

    /**
     * add the nontransmitted genotype
     * 
     * @param genotype
     */
    public void setFamilyID(String fid) {
        this.familyID = fid;
    }

    public void setPseudoPersonID(String pid) {
        this.pseudopersonID = pid;
    }

    public ArrayList<String> getGenotype() {
		ArrayList<String> sub = NewIt.newArrayList();
		for (int i = 0; i < numMarkers; i++) {
			StringBuffer sb = new StringBuffer();
			sb.append(alleles[0][i]);
			sb.append(alleles[1][i]);
			sub.add(sb.toString());
		}
		return sub; 
    }
 
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

    public String getFamilyID() {
        return familyID;
    }

    public String getPseudoPersonID() {
        return pseudopersonID;
    }

    public int getNumMarkers() {
        return numMarkers;
    }

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

    public ArrayList<String> getPseudoGenotype(int[] subsetMarker) {
		ArrayList<String> sub = NewIt.newArrayList();
		for (int i = 0; i < subsetMarker.length; i++) {
			StringBuffer sb = new StringBuffer();
			sb.append(alleles[0][subsetMarker[i]]);
			sb.append(alleles[1][subsetMarker[i]]);
			sub.add(sb.toString());
		}
		return sub; 
    }
}
