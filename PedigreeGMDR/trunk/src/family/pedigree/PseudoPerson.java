package family.pedigree;

import java.util.ArrayList;

/**
 * stores the genotypes of each individual, and, furthermore, store the nontransmitted genotype of this individual. this
 * class is not thread safe (untested)
 * 
 * @author Guo-Bo Chen
 */
public class PseudoPerson {

    private String familyID;
    private String pseudopersonID;
    private String momID;
    private String dadID;
    private int gender;
    private int affectedStatus;
    private ArrayList pseudoGenotype;
    private int numMarkers;

    PseudoPerson(int numMarkers) {
        this.numMarkers = numMarkers;
        pseudoGenotype = new ArrayList();
    }

    public void pseudoGenotypeClear() {
        pseudoGenotype.clear();
    }

    public void addMarker(String genotype) {
        pseudoGenotype.add(new String(genotype));
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

    public ArrayList getGenotype() {
        return pseudoGenotype;
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

    public ArrayList getPseudoGenotype() {
        return pseudoGenotype;
    }
    
    public ArrayList getPseudoGenotype(int[] subsetMarker) {
    	if(subsetMarker.length == pseudoGenotype.size()) {
    		return pseudoGenotype;
    	} else {
    		ArrayList sub = new ArrayList();
    		for(int i = 0; i < subsetMarker.length; i++) {
    			sub.add((String) pseudoGenotype.get(subsetMarker[i]));
    		}
    		return sub;
    	}
    }
}
