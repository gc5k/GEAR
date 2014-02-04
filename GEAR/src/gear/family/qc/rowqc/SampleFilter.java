package gear.family.qc.rowqc;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;

import gear.CmdArgs;
import gear.ConstValues;
import gear.data.Person;
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

	protected int qualified_Unrelated;
	protected int qualified_Sib;
	protected int[] numSib;
	protected double[] status;
	protected double[] score;
	// protected double[] permuted_score;

	protected ArrayList<PersonIndex> PersonTable;// The indexing file records
	private ArrayList<Hukou> HukouBook;
	protected int[][] num_qualified;//
	protected boolean[][] keep;

	protected boolean[][] filter;

	private String[][] indKeep;
	private String[][] indExclude;

	public SampleFilter(PedigreeFile ped, MapFile map)
	{
		PedData = ped;
		MapData = map;
		PersonTable = NewIt.newArrayList();

		QC();

	}

	public void QC()
	{
		qualification();
		if (PersonTable.size() == 0)
		{
			Logger.printUserError("No individual is left for analysis.");
			System.exit(1);
		}
	}

	private void qualification()
	{
		ArrayList<Hukou> hukoubook = PedData.getHukouBook();
		HukouBook = NewIt.newArrayList();
		UniqueRecordSet<Family> families = PedData.getFamilies();
		num_qualified = new int[families.size()][2];
		filter = new boolean[families.size()][];

		for (Iterator<Hukou> e = hukoubook.iterator(); e.hasNext();)
		{
			Hukou hukou = e.next();
			Family family = families.get(hukou.getFamilyID());
			Person per = family.getPerson(hukou.getIndividualID());
			boolean hf = hardFilter(per);
			hukou.setAvailable(hf);
			if (!hf)
			{
				continue;
			}
			boolean isFounder = family.hasAncestor(per);
			HukouBook.add(hukou);
			PersonTable.add(new PersonIndex(per.getFamilyID(), per
					.getPersonID(), per, false, isFounder));
		}

		num_qualified = new int[families.size()][2];
		filter = new boolean[families.size()][];
		int c = 0;
		for (int familyIdx = 0; familyIdx < families.size(); ++familyIdx)
		{
			Family family = families.get(familyIdx);
			filter[c] = new boolean[family.size()];

			// filter_family

			int cc = 0;
			for (int personIdx = 0; personIdx < family.size(); ++personIdx)
			{
				Person per = family.getPerson(personIdx);
				boolean hf = hardFilter(per);
				if (!hf)
				{
					filter[c][cc++] = hf;
					continue;
				}

				filter[c][cc++] = hf;
				boolean isFounder = family.hasAncestor(per);
				if (isFounder)
				{
					num_qualified[c][1]++;
				} else
				{
					num_qualified[c][0]++;
				}
			}
			c++;
		}

		for (int i = 0; i < num_qualified.length; i++)
		{
			qualified_Unrelated += num_qualified[i][0];
			qualified_Sib += num_qualified[i][1];
		}
	}

	protected boolean hardFilter(Person p)
	{
		boolean flag = true;

		if (CmdArgs.INSTANCE.keepFlag)
		{
			Logger.printUserLog("Reading kept individuals from '"
					+ CmdArgs.INSTANCE.keepFile + "'.");
			readKeepFile();
			flag = false;
			String fi = p.getFamilyID();
			String pi = p.getPersonID();
			for (int i = 0; i < indKeep[0].length; i++)
			{
				if (fi.compareTo(indKeep[0][i]) == 0
						&& pi.compareTo(indKeep[1][i]) == 0)
				{
					flag = true;
					break;
				}
			}
		} else if (CmdArgs.INSTANCE.removeFlag)
		{
			Logger.printUserLog("Reading removed individuals from '"
					+ CmdArgs.INSTANCE.removeFile + "'.");
			readRemoveFile();
			String fi = p.getFamilyID();
			String pi = p.getPersonID();
			for (int i = 0; i < indExclude[0].length; i++)
			{
				if (fi.compareTo(indExclude[0][i]) == 0
						&& pi.compareTo(indExclude[1][i]) == 0)
				{
					flag = false;
					break;
				}
			}
		}

		if (flag && CmdArgs.INSTANCE.keep_maleFlag)
		{
			return flag = p.getGender() == 1 ? true : false;
		}
		if (flag && CmdArgs.INSTANCE.keep_femaleFlag)
		{
			return flag = p.getGender() == 2 ? true : false;
		}
		if (flag && CmdArgs.INSTANCE.ex_nosexFlag)
		{
			return flag = (p.getGender() == 1 || p.getGender() == 2) ? true
					: false;
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
		BufferedReader reader = FileUtil
				.FileOpen(CmdArgs.INSTANCE.keepFile);
		String line = null;
		ArrayList<String> famList = NewIt.newArrayList();
		ArrayList<String> indList = NewIt.newArrayList();

		try
		{
			while ((line = reader.readLine()) != null)
			{
				String[] l = line.split(ConstValues.WHITESPACE_DELIMITER);
				if (l.length < 2)
					continue;
				famList.add(l[0]);
				indList.add(l[1]);
			}
		} catch (IOException e)
		{
			Logger.handleException(e,
					"An exception occurred when reading the kept-individual file '"
							+ CmdArgs.INSTANCE.keepFile + "'.");
		}
		indKeep = new String[2][];
		indKeep[0] = (String[]) famList.toArray(new String[0]);
		indKeep[1] = (String[]) indList.toArray(new String[0]);
	}

	private void readRemoveFile()
	{
		BufferedReader reader = FileUtil
				.FileOpen(CmdArgs.INSTANCE.removeFile);
		ArrayList<String> famList = NewIt.newArrayList();
		ArrayList<String> indList = NewIt.newArrayList();
		String line = null;
		try
		{
			while ((line = reader.readLine()) != null)
			{
				String[] l = line.split(ConstValues.WHITESPACE_DELIMITER);
				if (l.length < 2)
					continue;
				famList.add(l[0]);
				indList.add(l[1]);
			}
		} catch (IOException e)
		{
			Logger.handleException(e,
					"An exception occurred when reading the removed-individual file '"
							+ CmdArgs.INSTANCE.removeFile + "'.");
		}
		indExclude = new String[2][];
		indExclude[0] = (String[]) famList.toArray(new String[0]);
		indExclude[1] = (String[]) indList.toArray(new String[0]);

	}
}
