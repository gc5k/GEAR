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
public class InputDataSet2 implements SubjectOrder
{
	public InputDataSet2()
	{
		
	}

	public void addFile(String subFile)
	{
		PhenotypeFile subjectIDFile = new PhenotypeFile(subFile);
		Logger.printUserLog("Read " + subjectIDFile.getNumberOfSubjects() + " samples from '" + subjectIDFile.getFileName() + "'.");
		fileList.add(subjectIDFile);
		fileNameList.add(subFile);
		makeSubjectOrderConsistent1(subjectIDFile);
	}

	public void addFile(String subFile, boolean flag)
	{
		PhenotypeFile subjectIDFile = new PhenotypeFile(subFile, ConstValues.NO_HEADER);
		Logger.printUserLog("Read " + subjectIDFile.getNumberOfSubjects() + " samples from '" + subjectIDFile.getFileName() + "'.");
		fileList.add(subjectIDFile);
		fileNameList.add(subFile);
		int[] tIdx = new int[subjectIDFile.getNumberOfTraits()];
		for (int i = 0; i < tIdx.length; i++) tIdx[i] = i;
		makeSubjectOrderConsistent1(subjectIDFile, tIdx);
	}
	
	public void addFile(String subFile, int[] tIdx)
	{
		PhenotypeFile subjectIDFile = new PhenotypeFile(subFile, ConstValues.NO_HEADER);
		Logger.printUserLog("Read " + subjectIDFile.getNumberOfSubjects() + " samples from '" + subjectIDFile.getFileName() + "'.");
		fileList.add(subjectIDFile);
		fileNameList.add(subFile);
		makeSubjectOrderConsistent1(subjectIDFile, tIdx);
	}

	public void LineUpFiles()
	{
		matchedIDAcrossFiles();
		for (int i = 0; i < fileList.size(); i++)
		{
			PhenotypeFile pf = fileList.get(i);
			int[] idx = new int[effectiveSubList.size()];
			ArrayList<SubjectID> SID = NewIt.newArrayList();
			for (int j = 0; j < effectiveSubList.size(); j++)
			{
				SubjectID sid = effectiveSubList.get(j);
				int subidx = pf.getSubjectIndex(sid);
				idx[j] = subidx;
				SID.add(sid);
			}
			sampleIdx.add(idx);
			sampleID.add(SID);
		}

		if (effectiveSubList.size() == 0) {
			Logger.printUserLog("No individuals were remained for analysis. GEAR quit.");
			System.exit(1);
		}
		if (fileList.size() > 0) {
			Logger.printUserLog(effectiveSubList.size() + " samples were matched across " + fileList.size() + " files.");			
		}
		Logger.printUserLog("");
	}

	private void matchedIDAcrossFiles()
	{
		ArrayList<SubjectID> sList = NewIt.newArrayList();
		Set<SubjectID> keySet = id2Idx.keySet();
		for(Iterator<SubjectID> e = keySet.iterator(); e.hasNext(); )
		{
			SubjectID sid = e.next();
			if (id2Idx.get(sid) == fileList.size())
			{
				sList.add(sid);
			}
		}
		if (sList.size() < 1)
		{
			Logger.printUserLog("No samples were lined up!");
			Logger.printUserLog("GEAR quitted.\n");
			System.exit(1);
		}
		
		PhenotypeFile pf = fileList.get(0);
		
		for (int i = 0; i < pf.getNumberOfSubjects(); i++)
		{
			SubjectID sid = pf.getSubjectID(i);
			if (sList.contains(sid)) effectiveSubList.add(sid);
		}	
	}

	@Override
	public int getNumberOfSubjects()
	{
		return effectiveSubList.size();
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

	public ArrayList<SubjectID> getMatchedSubjectID(int fileIdx)
	{
		return sampleID.get(fileIdx);
	}

	@Override
	public void swapSubjects(int subjectIdx1, int subjectIdx2)
	{
	}

	public double[] getVariable(int fileIndex, int subjectIdx, int[] covIdx)
	{
		PhenotypeFile pf = fileList.get(fileIndex);
		double[] vars = new double[covIdx.length];
		for (int i = 0; i < covIdx.length; i++)
		{
			vars[i] = (double) pf.getPhenotype(subjectIdx, covIdx[i]);
		}
		return vars;
	}

	public double getVariable(int fileIndex, int subjectIdx, int covIdx)
	{
		PhenotypeFile pf = fileList.get(fileIndex);
		return (double) pf.getPhenotype(subjectIdx, covIdx);
	}

	public boolean isVariableMissing(int fileIndex, int subjectIdx, int[] covIdx)
	{
		boolean flag = true;
		PhenotypeFile pf = fileList.get(fileIndex);
		for (int i = 0; i < covIdx.length; i++)
		{
			flag = pf.isMissing(subjectIdx, covIdx[i]);
		}
		return flag;
	}

	public boolean isVariableMissing(int fileIndex, int subjectIdx, int covIdx)
	{
		PhenotypeFile pf = fileList.get(fileIndex);
		return pf.isMissing(subjectIdx, covIdx);
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

	private void makeSubjectOrderConsistent1(SubjectOrder order, int[] tIdx)
	{
		PhenotypeFile pf = fileList.get(fileList.size()-1);
		for (int subjectIdx1 = 0; subjectIdx1 < order.getNumberOfSubjects() ; ++subjectIdx1)
		{
			SubjectID subjectID = order.getSubjectID(subjectIdx1);
			boolean flag = false;
			for (int j = 0; j < tIdx.length; j++)
			{
				flag = pf.isMissing(pf.getSubjectIndex(subjectID), tIdx[j]);
			}
			if (flag) continue;

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

	public int getNumberOfTraits()
	{
		return fileList.get(1).getNumberOfTraits();
	}

	public int getFileSampleSize(int idx)
	{
		PhenotypeFile pf = fileList.get(idx);
		return pf.getNumberOfSubjects();
	}

	public int[] getMatchedSubjectIdx(int idx)
	{
		int[] sid = sampleIdx.get(idx);
		return sid;
	}

	public ArrayList<SubjectID> getMatchSubjetList()
	{
		return effectiveSubList;
	}

	private ArrayList<SubjectID> effectiveSubList = NewIt.newArrayList();
	private SubjectOrder sharedSubjectOrder;

	private HashMap<SubjectID, Integer> id2Idx = new HashMap<SubjectID, Integer>();
	private ArrayList<PhenotypeFile> fileList = NewIt.newArrayList();
	private ArrayList<int[]> sampleIdx = NewIt.newArrayList();
	private ArrayList<ArrayList<SubjectID>> sampleID = NewIt.newArrayList();
	private ArrayList<String> fileNameList = NewIt.newArrayList();
}
