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
import gear.family.pedigree.file.PedigreeFile;
import gear.family.plink.PLINKParser;
import gear.subcommands.CommandArguments;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;

/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */

public class SampleFilter {
	protected PedigreeFile PedData;
	protected PLINKParser plinkParser;

	protected String kFile;
	protected String rFile;
	protected int filterType = 0; // 0 for no filter, 1 for keep, 2 for remove

	protected double[] status;
	protected double[] score;

	protected ArrayList<PersonIndex> PersonTable;// The indexing file records
	private ArrayList<Hukou> HukouBook;
	protected int[][] num_qualified;//
	protected boolean[][] keep;

	protected boolean[][] filter;

	private HashSet<SubjectID> subSet;

	public SampleFilter(PedigreeFile ped) {
		PedData = ped;
		PersonTable = NewIt.newArrayList();

		qualification();

	}

	public SampleFilter(PedigreeFile ped, CommandArguments cmdArgs) {
		PedData = ped;
		PersonTable = NewIt.newArrayList();

		if (cmdArgs.isKeepFile()) {
			kFile = cmdArgs.getKeepFile();
			readKeepFile();
			filterType = 1;
		} else if (cmdArgs.isRemoveFile()) {
			rFile = cmdArgs.getRemoveFile();
			readRemoveFile();
			filterType = 2;
		}
		
		qualification();
	}

	public SampleFilter(PedigreeFile ped, ArrayList<SubjectID> sID) {
		PedData = ped;
		PersonTable = NewIt.newArrayList();
		filterType = 1;
		subSet = NewIt.newHashSet();
		for (SubjectID sid : sID) {
			subSet.add(sid);
		}
		qualification();
	}

	private void qualification() {
		ArrayList<Hukou> hukoubook = PedData.getHukouBook();
		HukouBook = NewIt.newArrayList();
		UniqueRecordSet<Family> families = PedData.getFamilies();

		for (Iterator<Hukou> e = hukoubook.iterator(); e.hasNext();) {
			Hukou hukou = e.next();
			Family family = families.get(hukou.getFamilyID());
			Person per = family.getPerson(hukou.getIndividualID());
			boolean isKeep = keep(per);
			hukou.setAvailable(isKeep);
			if (!isKeep) {
				continue;
			}
			boolean isFounder = family.hasAncestor(per);
			HukouBook.add(hukou);
			PersonTable.add(new PersonIndex(per, false, isFounder));
		}
		if (PersonTable.size() == 0) {
			Logger.printUserLog("No individuals were remained for analysis. GEAR quit.");
			System.exit(1);
		}

		Logger.printUserLog(PersonTable.size() + " individuals were matched for analysis.");
	}

	protected boolean keep(Person p) {
		boolean flag = true;

		if (filterType == 0) {
			return flag;
		} else if (filterType == 1) {
			flag = subSet.contains(new SubjectID(p.getFamilyID(), p.getPersonID()));
		} else if (filterType == 2) {
			flag = !subSet.contains(new SubjectID(p.getFamilyID(), p.getPersonID()));
		}
		return flag;
	}

	public int SampleSize() {
		return PersonTable.size();
	}

	public ArrayList<PersonIndex> getSample() {
		return PersonTable;
	}

	public ArrayList<Hukou> getHukouBook() {
		return HukouBook;
	}

	private void readKeepFile() {
		BufferedReader reader = FileUtil.FileOpen(kFile);
		String line = null;
		subSet = NewIt.newHashSet();

		try {
			while ((line = reader.readLine()) != null) {
				String[] l = line.split(ConstValues.WHITESPACE_DELIMITER);
				if (l.length < 2)
					continue;
				subSet.add(new SubjectID(l[0], l[1]));
			}
		} catch (IOException e) {
			Logger.handleException(e, "An exception occurred when reading the keep file '" + kFile + "'.");
		}
		Logger.printUserLog("Reading " + subSet.size() + " individuals from keep-individual file '" + kFile + "'.");
	}

	private void readRemoveFile() {
		BufferedReader reader = FileUtil.FileOpen(rFile);
		String line = null;
		subSet = NewIt.newHashSet();

		try {
			while ((line = reader.readLine()) != null) {
				String[] l = line.split(ConstValues.WHITESPACE_DELIMITER);
				if (l.length < 2)
					continue;
				subSet.add(new SubjectID(l[0], l[1]));
			}
		} catch (IOException e) {
			Logger.handleException(e,
					"An exception occurred when reading the removed-individual file '" + rFile + "'.");
		}
		Logger.printUserLog("Reading " + subSet.size() + " individuals from remove-individal file '" + rFile + "'.");
	}
}
