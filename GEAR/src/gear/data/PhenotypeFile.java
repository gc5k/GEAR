package gear.data;

import gear.ConstValues;
import gear.util.BufferedReader;
import gear.util.IDIndexMap;
import gear.util.Logger;

import java.util.ArrayList;

public class PhenotypeFile implements SubjectOrder
{
	public PhenotypeFile(String fileName)
	{
		this.fileName = fileName;
		
		BufferedReader reader = BufferedReader.openTextFile(fileName, "phenotype");

		String[] tokens = reader.readTokensAtLeast(3);
		
		if (tokens == null)
		{
			Logger.printUserError("The phenotype file '" + fileName + "' is empty.");
			System.exit(1);
		}
		
		int numCols = tokens.length;
		ArrayList<float[]> phenotypes = new ArrayList<float[]>();
		ArrayList<boolean[]> isMissing = new ArrayList<boolean[]>();

		do
		{
			SubjectID subjectID = new SubjectID(/*famID*/tokens[0], /*indID*/tokens[1]);

			if (!subjectMap.add(subjectID))
			{
				reader.errorPreviousLine("Subject " + subjectID + " is duplicated.");
			}
			
			float[] phenotypesOfThisSubject = new float[numCols - 2];
			boolean[] isMissingOfThisSubject = new boolean[numCols - 2];
			
			for( int traitIdx = 0; traitIdx < phenotypesOfThisSubject.length; traitIdx++)
			{
				String pheValStr = tokens[2 + traitIdx];
				if (ConstValues.isNA(pheValStr))
				{
					isMissingOfThisSubject[traitIdx] = true;
				}
				else
				{
					try
					{
						phenotypesOfThisSubject[traitIdx] = Float.parseFloat(pheValStr);
					}
					catch (NumberFormatException e)
					{
						reader.errorPreviousLine("'" + pheValStr + "' is not a valid phenotype value. It should be a floating point number.");
					}
				}						
			}
			
			phenotypes.add(phenotypesOfThisSubject);
			isMissing.add(isMissingOfThisSubject);
		} while ((tokens = reader.readTokens(numCols)) != null);
		
		reader.close();

		this.phenotypes = phenotypes.toArray(new float[0][]);
		this.isMissing = isMissing.toArray(new boolean[0][]);
	}
	
	public String getFileName()
	{
		return fileName;
	}
	
	@Override
	public String getSubjectOrderName()
	{
		return "file '" + getFileName() + "'";
	}
	
	@Override
	public int getNumberOfSubjects()
	{
		return phenotypes.length;
	}
	
	@Override
	public int getSubjectIndex(SubjectID subjectID)
	{
		return subjectMap.getIndex(subjectID);
	}
	
	@Override
	public SubjectID getSubjectID(int subjectIdx)
	{
		return subjectMap.getID(subjectIdx);
	}
	
	@Override
	public void swapSubjects(int subjectIdx1, int subjectIdx2)
	{
		subjectMap.swapEntries(subjectIdx1, subjectIdx2);
		
		float[] tmpPhenotypes = phenotypes[subjectIdx1];
		phenotypes[subjectIdx1] = phenotypes[subjectIdx2];
		phenotypes[subjectIdx2] = tmpPhenotypes;
		
		boolean[] tmpIsMissing = isMissing[subjectIdx1];
		isMissing[subjectIdx1] = isMissing[subjectIdx2];
		isMissing[subjectIdx2] = tmpIsMissing;
	}
	
	public int getNumberOfTraits()
	{
		return phenotypes[0].length;
	}
	
	public boolean isMissing(int subjectIdx, int traitIdx)
	{
		return isMissing[subjectIdx][traitIdx];
	}
	
	public float getPhenotype(int subjectIdx, int traitIdx)
	{
		if (isMissing(subjectIdx, traitIdx))
		{
			Logger.warnInternalBug("Retrieving a missing phenotype value.");
		}
		return phenotypes[subjectIdx][traitIdx];
	}
	
	private String fileName;
	private IDIndexMap<SubjectID> subjectMap = new IDIndexMap<SubjectID>();
	private float[][] phenotypes;
	private boolean[][] isMissing;
}
