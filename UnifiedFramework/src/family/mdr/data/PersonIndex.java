package family.mdr.data;

import family.mdr.arsenal.MDRConstant;
import family.pedigree.genotype.BPerson;

public class PersonIndex {

	private String FamilyID;
	private String IndividualID;
	private BPerson person;
	private int group;
	private double permutatedScore;
	
	public PersonIndex(String fid, String id, BPerson p) {
		FamilyID = new String(fid);
		IndividualID = id;
		person = p;
	}

	public String getFamilyID() {
		return FamilyID;
	}

	public String getIndividualID() {
		return IndividualID;
	}
	
	public BPerson getPerson() {
		return person;
	}
	
	public double getScore() {
		return permutatedScore;
	}
	
	public void setPermutedScore(double s) {
		permutatedScore = s;
	}
	
	public void setGroup(int g) {
		group = g;
	}
	
	public int getGroup() {
		return group;
	}
	
	public String getGenotype(int[] idx) {
		StringBuffer sb = new StringBuffer();
		for(int i = 0; i < idx.length; i++) {
			sb.append(person.getGenotypeScoreString(idx[i]));
			sb.append(MDRConstant.seperator);
		}
		return sb.toString();
	}
}
