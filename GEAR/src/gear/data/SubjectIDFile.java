package gear.data;

import gear.util.BufferedReader;
import gear.util.IDIndexMap;

public class SubjectIDFile implements SubjectOrder
{
	public SubjectIDFile(String fileName)
	{
		this.fileName = fileName;
		BufferedReader reader = BufferedReader.openTextFile(fileName, "subject-ID");
		String[] tokens;
		while ((tokens = reader.readTokensAtLeast(2)) != null)
		{
			SubjectID id = new SubjectID(/*famID*/tokens[0], /*indID*/tokens[1]);
			if (!subjectMap.add(id))
			{
				reader.errorPreviousLine("Subject " + id + " is duplicated.");
			}
		}
		
		if (subjectMap.getNumberOfEntries() == 0)
		{
			reader.errorPreviousLine("The file is empty.");
		}
		
		reader.close();
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
		return subjectMap.getNumberOfEntries();
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
	}
	
	private String fileName;
	private IDIndexMap<SubjectID> subjectMap = new IDIndexMap<SubjectID>();
}
