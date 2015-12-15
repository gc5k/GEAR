package gear.data;

import gear.ConstValues;
import gear.util.Logger;

/**
 * This class acts as a builder and facade (see design patterns [Gang of Four]
 * for reference) to provide a simplified interface for accessing common data
 * files.
 * 
 * @author Zhixiang
 *
 */
public class InputDataSet implements SubjectOrder
{
	public void readSubjectIDFile(String fileName)
	{
		subjectIDFile = new SubjectIDFile(fileName);
		makeSubjectOrderConsistent(subjectIDFile);
	}

	public void readPhenotypeFile(String fileName)
	{
		phenotypeFile = new PhenotypeFile(fileName, ConstValues.NO_HEADER);
		makeSubjectOrderConsistent(phenotypeFile);
	}
	
	@Override
	public int getNumberOfSubjects()
	{
		return sharedSubjectOrder.getNumberOfSubjects();
	}
	
	@Override
	public String getSubjectOrderName()
	{
		return "final data used for computation";
	}

	@Override
	public int getSubjectIndex(SubjectID subjectID)
	{
		return sharedSubjectOrder.getSubjectIndex(subjectID);
	}

	@Override
	public SubjectID getSubjectID(int subjectIdx)
	{
		return sharedSubjectOrder.getSubjectID(subjectIdx);
	}

	@Override
	public void swapSubjects(int subjectIdx1, int subjectIdx2)
	{
		if (subjectIDFile != null)
		{
			subjectIDFile.swapSubjects(subjectIdx1, subjectIdx2);
		}
		
		if (phenotypeFile != null)
		{
			phenotypeFile.swapSubjects(subjectIdx1, subjectIdx2);
		}
	}
	
	public int getNumberOfTraits()
	{
		return phenotypeFile.getNumberOfTraits();
	}
	
	public float getPhenotype(int subjectIdx, int traitIdx)
	{
		return phenotypeFile.getPhenotype(subjectIdx, traitIdx);
	}
	
	public boolean isPhenotypeMissing(int subjectIdx, int traitIdx)
	{
		return phenotypeFile.isMissing(subjectIdx, traitIdx);
	}
	
	private void makeSubjectOrderConsistent(SubjectOrder order)
	{
		if (sharedSubjectOrder == null)
		{
			sharedSubjectOrder = order;
		}
		else if (sharedSubjectOrder != order)
		{
			makeSubjectOrdersConsistent(sharedSubjectOrder, order);
		}
	}
	
	private static void makeSubjectOrdersConsistent(SubjectOrder order1, SubjectOrder order2)
	{
		if (order1.getNumberOfSubjects() != order2.getNumberOfSubjects())
		{
			String msg = "";
			msg += "The number of subjects in " + order1.getSubjectOrderName() + " is " + order1.getNumberOfSubjects() + ", ";
			msg += "but the number of subjects in " + order2.getSubjectOrderName() + " is " + order2.getNumberOfSubjects() + ". ";
			msg += "They are inconsistent.";
			Logger.printUserError(msg);
			System.exit(1);
		}
		
		for (int subjectIdx1 = 0; subjectIdx1 < order1.getNumberOfSubjects() - 1; ++subjectIdx1)
		{
			SubjectID subjectID = order1.getSubjectID(subjectIdx1);
			int subjectIdx2 = order2.getSubjectIndex(subjectID);
			
			if (subjectIdx2 < 0 || subjectIdx2 >= order2.getNumberOfSubjects())
			{
				String msg = "";
				msg += "Subject " + subjectID + " appears in " + order1.getSubjectOrderName();
				msg += " but not in " + order2.getSubjectOrderName() + ".";
				Logger.printUserError(msg);
				System.exit(1);
			}
			
			if (subjectIdx1 != subjectIdx2)
			{
				order2.swapSubjects(subjectIdx1, subjectIdx2);
			}
		}
	}
	
	private SubjectOrder sharedSubjectOrder;
	private SubjectIDFile subjectIDFile;
	private PhenotypeFile phenotypeFile;
}
