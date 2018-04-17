package gear.family.qc.rowqc;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;

import gear.ConstValues;
import gear.data.Person;
import gear.data.SubjectID;
import gear.data.Family;
import gear.data.UniqueRecordSet;
import gear.family.pedigree.Hukou;
import gear.family.pedigree.PersonIndex;
import gear.family.pedigree.file.MapFile;
import gear.family.pedigree.file.PedigreeFile;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;

/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */

public class SampleFilter
{

	protected MapFile MapData;
	protected PedigreeFile PedData;

	protected String kFile;
	protected String rFile;
	protected int filterType = 0; //0 for no filter, 1 for keep, 2 for remove

	protected double[] status;
	protected double[] score;

	protected ArrayList<PersonIndex> PersonTable;// The indexing file records
	private ArrayList<Hukou> HukouBook;
	protected int[][] num_qualified;//
	protected boolean[][] keep;

	protected boolean[][] filter;

	private HashSet<SubjectID> subSet;

	public SampleFilter(PedigreeFile ped, MapFile map)
	{
		PedData = ped;
		MapData = map;
		PersonTable = NewIt.newArrayList();

		QC();

	}

	public SampleFilter(PedigreeFile ped, MapFile map, String keepFile, String removeFile)
	{
		PedData = ped;
		MapData = map;
		PersonTable = NewIt.newArrayList();

		if (keepFile != null)
		{
			kFile = keepFile;
			readKeepFile();
			filterType = 1;
		} else if (removeFile != null) {
			rFile = removeFile;
			readRemoveFile();
			filterType = 2;
		}
		QC();
	}

	public void QC()
	{
		qualification();
		if (PersonTable.size() == 0)
		{
			Logger.printUserError("No individuals are left for analysis.");
			System.exit(1);
		}
	}

	private void qualification()
	{
		ArrayList<Hukou> hukoubook = PedData.getHukouBook();
		HukouBook = NewIt.newArrayList();
		UniqueRecordSet<Family> families = PedData.getFamilies();

		for (Iterator<Hukou> e = hukoubook.iterator(); e.hasNext();)
		{
			Hukou hukou = e.next();
			Family family = families.get(hukou.getFamilyID());
			Person per = family.getPerson(hukou.getIndividualID());
			boolean isKeep = keep(per);
			hukou.setAvailable(isKeep);
			if (!isKeep)
			{
				continue;
			}
			boolean isFounder = family.hasAncestor(per);
			HukouBook.add(hukou);
			PersonTable.add(new PersonIndex(per.getFamilyID(), per
					.getPersonID(), per, false, isFounder));
		}
		Logger.printUserLog(PersonTable.size() + " individuals were remained for analysis.");
	}

	protected boolean keep(Person p)
	{
		boolean flag = true;

		if (filterType == 0)
		{
			return flag;
		}
		else if (filterType == 1)
		{
			flag = subSet.contains(new SubjectID(p.getFamilyID(), p.getPersonID()));
		}
		else if (filterType == 2)
		{
			flag = !subSet.contains(new SubjectID(p.getFamilyID(), p.getPersonID()));
		}
		return flag;
	}

	public int getNumberMarker()
	{
		return MapData.getMarkerNumber();
	}

	public MapFile getMapFile()
	{
		return MapData;
	}

	public int SampleSize()
	{
		return PersonTable.size();
	}

	public ArrayList<PersonIndex> getSample()
	{
		return PersonTable;
	}

	public ArrayList<Hukou> getHukouBook()
	{
		return HukouBook;
	}

	private void readKeepFile()
	{
		BufferedReader reader = FileUtil.FileOpen(kFile);
		String line = null;
		subSet = NewIt.newHashSet();

		try
		{
			while ((line = reader.readLine()) != null)
			{
				String[] l = line.split(ConstValues.WHITESPACE_DELIMITER);
				if (l.length < 2)
					continue;
				subSet.add(new SubjectID(l[0], l[1]));
			}
		} catch (IOException e)
		{
			Logger.handleException(e,
					"An exception occurred when reading the keep file '"
							+ kFile + "'.");
		}
		Logger.printUserLog("Reading " + subSet.size() + " individuals from keep-individual file '"+ kFile + "'.");
	}

	private void readRemoveFile()
	{
		BufferedReader reader = FileUtil.FileOpen(rFile);
		String line = null;
		subSet = NewIt.newHashSet();

		try
		{
			while ((line = reader.readLine()) != null)
			{
				String[] l = line.split(ConstValues.WHITESPACE_DELIMITER);
				if (l.length < 2)
					continue;
				subSet.add(new SubjectID(l[0], l[1]));
			}
		} catch (IOException e)
		{
			Logger.handleException(e,
					"An exception occurred when reading the removed-individual file '"
							+ rFile + "'.");
		}
		Logger.printUserLog("Reading " + subSet.size() + " individuals from remove-individal file '"+ rFile + "'.");
	}
}
