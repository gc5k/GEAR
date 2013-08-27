package gear.family.pedigree.phenotype;

import gear.util.NewIt;

import java.util.ArrayList;

/**
 * stores the phenotypic outcomes of each individual. this class is not thread
 * safe (untested)
 * 
 * @author Guo-Bo Chen
 */
public class Subject
{

	private String familyID;
	private String subjectID;
	private ArrayList<String> traits;
	private int numTraits;

	public Subject(int num)
	{
		traits = NewIt.newArrayList();
		numTraits = num;
	}

	public void AddTrait(final String t)
	{
		traits.add(t);
	}

	public int getNumberofTraits()
	{
		return numTraits;
	}

	public ArrayList<String> getTraits()
	{
		return traits;
	}

	public void setFamilyID(String FID)
	{
		familyID = FID;
	}

	public String getFamilyID()
	{
		return familyID;
	}

	public void setSubjectID(String ID)
	{
		subjectID = ID;
	}

	public String getSubjectID()
	{
		return subjectID;
	}
}
