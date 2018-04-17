package gear.family.pedigree;

import gear.data.Person;
import gear.data.SubjectID;

public class PersonIndex
{
	private SubjectID subID;
	private Person person;
	private boolean pseudo;
	private boolean isFounder;
	private int group;
	private double permutatedScore;

	public PersonIndex(String FamilyID, String IndividualID, Person person,
			boolean pseudo, boolean isFounder)
	{
		subID = new SubjectID(FamilyID, IndividualID);
		this.person = person;

		this.pseudo = pseudo;
		this.isFounder = isFounder;
	}

	public String getFamilyID()
	{
		return subID.getFamilyID();
	}

	public String getIndividualID()
	{
		return subID.getIndividualID();
	}

	public Person getPerson()
	{
		return person;
	}

	public double getScore()
	{
		return permutatedScore;
	}

	public void setPermutedScore(double s)
	{
		permutatedScore = s;
	}

	public void setGroup(int g)
	{
		group = g;
	}

	public int getGroup()
	{
		return group;
	}

	public boolean isPseudo()
	{
		return pseudo;
	}

	public boolean isFounder()
	{
		return isFounder;
	}

	public String getGenotype(int[] idx)
	{
		StringBuffer sb = new StringBuffer();
		for (int i = 0; i < idx.length; i++)
		{
			sb.append(person.getGenotypeScoreString(idx[i]));
			sb.append(",");
		}
		return sb.toString();
	}
}
