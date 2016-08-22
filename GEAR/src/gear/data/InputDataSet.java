package gear.data;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;

import gear.ConstValues;
import gear.util.Logger;
import gear.util.NewIt;

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
	public InputDataSet()
	{
		
	}

	public InputDataSet(String subFile, String pheFile, int phenoIdx)
	{
		subjectIDFile = new SubjectIDFile(subFile);
		Logger.printUserLog("Read " + subjectIDFile.getNumberOfSubjects() + " samples from '" + subjectIDFile.getFileName() + "'.");
		makeSubjectOrderConsistent1(subjectIDFile);

		phenotypeFile = new PhenotypeFile(pheFile, ConstValues.NO_HEADER);
		Logger.printUserLog("Read " + phenotypeFile.getNumberOfSubjects() + " samples for " + phenotypeFile.getNumberOfTraits() + " traits from '" + phenotypeFile.getFileName() + "'.");
		if ((phenoIdx + 1) > phenotypeFile.getNumberOfTraits())
		{
			Logger.printUserLog((phenoIdx + 1) + " is incorrect phenotype index.");
			Logger.printUserLog("GEAR quitted.");
			System.exit(1);
		}

		makeSubjectOrderConsistent1(phenotypeFile);

		fileSize = 2;
		getMatchedID();
		lineupPhe(phenoIdx);
	}

	public InputDataSet(String subFile, String pheFile, String cFile, int phenoIdx, int[] covarIdx)
	{
		subjectIDFile = new SubjectIDFile(subFile);
		Logger.printUserLog("Read " + subjectIDFile.getNumberOfSubjects() + " samples from '" + subjectIDFile.getFileName() + "'.");
		makeSubjectOrderConsistent1(subjectIDFile);

		phenotypeFile = new PhenotypeFile(pheFile, ConstValues.NO_HEADER);
		Logger.printUserLog("Read " + phenotypeFile.getNumberOfSubjects() + " samples for " + phenotypeFile.getNumberOfTraits() + " trait(s) from '" + phenotypeFile.getFileName() + "'.");
		makeSubjectOrderConsistent1(phenotypeFile);

		covFile = new PhenotypeFile(cFile, ConstValues.NO_HEADER);
		Logger.printUserLog("Read " + covFile.getNumberOfSubjects() + " samples for " + covFile.getNumberOfTraits() + " covariate(s) from '" + covFile.getFileName() + "'.");
		makeSubjectOrderConsistent1(covFile);

		fileSize = 3;
		getMatchedID();
		lineupPheCov(phenoIdx, covarIdx);
	}

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

	public void readCovFile(String fileName)
	{
		covFile = new PhenotypeFile(fileName, ConstValues.NO_HEADER);
		makeSubjectOrderConsistent1(covFile);
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

	public float[] getCov(int subjectIdx, int[] covIdx)
	{
		float[] covs = new float[covIdx.length];
		for(int i = 0; i < covIdx.length; i++)
		{
			covs[i] = covFile.getPhenotype(subjectIdx, covIdx[i]);
		}
		return covs;
	}
	
	public float getCovariate(int subjectIdx, int cIdx)
	{
		return covFile.getPhenotype(subjectIdx, cIdx);
	}

	public boolean isCovariateMissing(int subjectIdx, int cIdx)
	{
		return covFile.isMissing(subjectIdx, cIdx);
	}
	
	public float getPhenotype(int subjectIdx, int traitIdx)
	{
		return phenotypeFile.getPhenotype(subjectIdx, traitIdx);
	}

	public boolean isPhenotypeMissing(int subjectIdx, int traitIdx)
	{
		return phenotypeFile.isMissing(subjectIdx, traitIdx);
	}

	private void makeSubjectOrderConsistent1(SubjectOrder order)
	{
		for (int subjectIdx1 = 0; subjectIdx1 < order.getNumberOfSubjects() ; ++subjectIdx1)
		{
			SubjectID subjectID = order.getSubjectID(subjectIdx1);
			if (id2Idx.containsKey(subjectID))
			{
				Integer cnt = id2Idx.get(subjectID);
				cnt++;
				id2Idx.put(subjectID, cnt);
			}
			else
			{
				id2Idx.put(subjectID, 1);
			}
		}
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

	private void getMatchedID()
	{
		Set<SubjectID> keySet = id2Idx.keySet();
		for(Iterator<SubjectID> e = keySet.iterator(); e.hasNext(); )
		{
			SubjectID sid = e.next();
			if (id2Idx.get(sid) == fileSize)
			{
				subList.add(sid);
			}
		}
		deepSubList = subList;
		if (deepSubList.size() < 5)
		{
			Logger.printUserLog("Only " + deepSubList.size() + " samples lined up. Too few!");
			Logger.printUserLog("GEAR quitted.\n");
			System.exit(1);
		}
	}

	public void lineupPheCov(int pheIdx, int[] covIdx)
	{
		ArrayList<SubjectID> tmpList = NewIt.newArrayList();
		for (int i = 0; i < subList.size(); i++)
		{
			SubjectID sid = subList.get(i);
			int subidx = phenotypeFile.getSubjectIndex(sid);
			int c_subidx = covFile.getSubjectIndex(sid);
			if (phenotypeFile.isMissing(subidx, pheIdx)) continue;

			for (int j = 0; j < covIdx.length; j++)
			{
				if (covFile.isMissing(c_subidx, covIdx[j])) continue;
			}
			tmpList.add(sid);
		}

		if (tmpList.size() < 5)
		{
			Logger.printUserLog("Only " + tmpList.size() + " lined up. Too few!");
			Logger.printUserLog("GEAR quitted.\n");
			System.exit(1);
		}
		Logger.printUserLog(tmpList.size() + " common samples were found across '" + subjectIDFile.getFileName() + "', '" +phenotypeFile.getFileName() +"' and '" + covFile.getFileName() + "'.");
		
///////line up according to the order in the subjectIDFile.
		deepSubList = NewIt.newArrayList();
		
		for(int i = 0; i < subjectIDFile.getNumberOfSubjects(); i++)
		{
			SubjectID sid = subjectIDFile.getSubjectID(i);
			if (tmpList.contains(sid)) deepSubList.add(sid);
		}
		Logger.printUserLog(deepSubList.size() + " samples were lined up as shown in '" + subjectIDFile.getFileName() + "'.");
		Logger.printUserLog("");

		SubSort();
		PhenoSubSort();
		CovSubSort();
	}

	public void lineupPhe(int pheIdx)
	{
		ArrayList<SubjectID> tmpList = NewIt.newArrayList();
		for (int i = 0; i < subList.size(); i++)
		{
			SubjectID sid = subList.get(i);
			int subidx = phenotypeFile.getSubjectIndex(sid);
			if ( phenotypeFile.isMissing(subidx, pheIdx) ) continue;
			tmpList.add(sid);
		}
		if (tmpList.size() < 5)
		{
			Logger.printUserLog("Only " + tmpList.size() + " lined up. Too few!");
			Logger.printUserLog("GEAR quitted\n");
			System.exit(1);
		}
		Logger.printUserLog(tmpList.size() + " " + "common samples were found between '" + subjectIDFile.getFileName() + "' and '" +phenotypeFile.getFileName() +"'.");

///////line up according to the order in the subjectIDFile.
		deepSubList = NewIt.newArrayList();
		
		for(int i = 0; i < subjectIDFile.getNumberOfSubjects(); i++)
		{
			SubjectID sid = subjectIDFile.getSubjectID(i);
			if (tmpList.contains(sid)) deepSubList.add(sid);
		}
		Logger.printUserLog(deepSubList.size() + " samples were lined up as shown in '" + subjectIDFile.getFileName() + "'.");
		Logger.printUserLog("");

		SubSort();
		PhenoSubSort();
	}

	private void SubSort()
	{
		ArrayList<Integer> subID = NewIt.newArrayList();
		for(int i = 0; i < deepSubList.size(); i++)
		{
			subID.add(subjectIDFile.getSubjectIndex(deepSubList.get(i)));
		}
		subIdx = new int[subID.size()];
		for(int i = 0; i < subID.size(); i++) subIdx[i] = subID.get(i);
	}

	private void PhenoSubSort()
	{
		ArrayList<Integer> pheID = NewIt.newArrayList();
		for(int i = 0; i < deepSubList.size(); i++)
		{
			pheID.add(phenotypeFile.getSubjectIndex(deepSubList.get(i)));
		}
		pheSubIdx = new int[pheID.size()];
		for(int i = 0; i < pheID.size(); i++) pheSubIdx[i] = pheID.get(i);
	}

	private void CovSubSort()
	{
		ArrayList<Integer> covID = NewIt.newArrayList();
		for(int i = 0; i < deepSubList.size(); i++)
		{
			covID.add(covFile.getSubjectIndex(deepSubList.get(i)));
		}
		covSubIdx = new int[covID.size()];
		for(int i = 0; i < covID.size(); i++) covSubIdx[i] = covID.get(i);
	}

	public int getSubjectFileSampleSize()
	{
		return subjectIDFile.getNumberOfSubjects();
	}

	public int getPhenotypeFileSampleSize()
	{
		return phenotypeFile.getNumberOfSubjects();
	}

	public int getCovFileSampleSize()
	{
		return covFile.getNumberOfSubjects();
	}

	public int getSampleSize()
	{
		return subList.size();
	}
	
	public int[] getMatchedSubIdx()
	{
		return subIdx;
	}
	
	public int[] getMatchedPheSubIdx()
	{
		return pheSubIdx;
	}
	
	public int[] getMatchedCovSubIdx()
	{
		return covSubIdx;
	}

	private int fileSize = 2;
	private ArrayList<SubjectID> subList = NewIt.newArrayList();
	private ArrayList<SubjectID> deepSubList = null;
	private SubjectOrder sharedSubjectOrder;

	private SubjectIDFile subjectIDFile;
	private int[] subIdx = null;

	private PhenotypeFile phenotypeFile;
	private int[] pheSubIdx = null;

	private PhenotypeFile covFile;
	private int[] covSubIdx = null;
	
	private HashMap<SubjectID, Integer> id2Idx = new HashMap<SubjectID, Integer>();
}
