package family.pedigree;

import family.pedigree.genotype.BPerson;

public class PersonIndex {

	private String FamilyID;
	private String IndividualID;
	private BPerson person;
	private boolean pseudo;
	private boolean isFounder;
	private int group;
	private double permutatedScore;
	
	public PersonIndex(String FamilyID, String IndividualID, BPerson person, boolean pseudo, boolean isFounder) {
		this.FamilyID = FamilyID;
		this.IndividualID = IndividualID;
		this.person = person;
		
		this.pseudo = pseudo;
		this.isFounder = isFounder;
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
	
	public boolean isPseudo() {
		return pseudo;
	}
	
	public boolean isFounder() {
		return isFounder;
	}

	public String getGenotype(int[] idx) {
		StringBuffer sb = new StringBuffer();
		for(int i = 0; i < idx.length; i++) {
			sb.append(person.getGenotypeScoreString(idx[i]));
			sb.append(",");
		}
		return sb.toString();
	}
}
