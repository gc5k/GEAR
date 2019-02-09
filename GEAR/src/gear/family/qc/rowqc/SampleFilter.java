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
	protected PedigreeFile pedData;
	protected PLINKParser plinkParser;
	
	public enum FilterType { NONE, KEEP_SAMPLES, REMOVE_SAMPLES, KEEP_FAMILIES, REMOVE_FAMILIES }
	protected FilterType filterType = FilterType.NONE;

	protected double[] status;
	protected double[] score;

	protected ArrayList<PersonIndex> personTable;// The indexing file records
	private ArrayList<Hukou> hukouBook;
	protected int[][] num_qualified;//
	protected boolean[][] keep;

	protected boolean[][] filter;

	private HashSet<SubjectID> subSet;
	private HashSet<String> famSet;

	public SampleFilter(PedigreeFile ped) {
		pedData = ped;
		personTable = NewIt.newArrayList();
		qualification();
	}
	
	public SampleFilter(CommandArguments cmdArgs) {
		personTable = NewIt.newArrayList();
		if (cmdArgs.isKeepFile()) {
			readIndividualFile(cmdArgs.getKeepFile(), "keep-individual");
			filterType = FilterType.KEEP_SAMPLES;
		} else if (cmdArgs.isRemoveFile()) {
			readIndividualFile(cmdArgs.getRemoveFile(), "remove-individual");
			filterType = FilterType.REMOVE_SAMPLES;
		} else if (cmdArgs.isKeepFamFile()) {
			readFamilyFile(cmdArgs.getKeepFamFile(), "keep-family");
			filterType = FilterType.KEEP_FAMILIES;
		} else if (cmdArgs.isRemoveFamFile()) {
			readFamilyFile(cmdArgs.getRemoveFamFile(), "remove-family");
			filterType = FilterType.REMOVE_FAMILIES;
		}
	}

	public SampleFilter(PedigreeFile ped, CommandArguments cmdArgs) {
		this(cmdArgs);
		pedData = ped;
		qualification();
	}

	public SampleFilter(PedigreeFile ped, ArrayList<SubjectID> sID) {
		pedData = ped;
		personTable = NewIt.newArrayList();
		filterType = FilterType.KEEP_SAMPLES;
		subSet = NewIt.newHashSet();
		for (SubjectID sid : sID) {
			subSet.add(sid);
		}
		qualification();
	}

	private void qualification() {
		ArrayList<Hukou> hukoubook = pedData.getHukouBook();
		hukouBook = NewIt.newArrayList();
		UniqueRecordSet<Family> families = pedData.getFamilies();

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
			hukouBook.add(hukou);
			personTable.add(new PersonIndex(per, false, isFounder));
		}
		if (personTable.size() == 0) {
			Logger.printUserLog("No individuals were remained for analysis. GEAR quit.");
			System.exit(1);
		}

		Logger.printUserLog(personTable.size() + " individuals were matched for analysis.");
	}

	protected boolean keep(Person p) {
		switch (filterType) {
		case NONE:
			return true;
		case KEEP_SAMPLES:
			return subSet.contains(new SubjectID(p.getFamilyID(), p.getPersonID()));
		case REMOVE_SAMPLES:
			return !subSet.contains(new SubjectID(p.getFamilyID(), p.getPersonID()));
		case KEEP_FAMILIES:
			return famSet.contains(p.getFamilyID());
		case REMOVE_FAMILIES:
			return !famSet.contains(p.getFamilyID());
		}
		return true;
	}

	public int SampleSize() {
		return personTable.size();
	}

	public ArrayList<PersonIndex> getSample() {
		return personTable;
	}

	public ArrayList<Hukou> getHukouBook() {
		return hukouBook;
	}

	private void readIndividualFile(String kFile, String opt) {
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
			Logger.handleException(e, "An exception occurred when reading the " + opt + " file '" + kFile + "'.");
		}
		Logger.printUserLog("Read " + subSet.size() + " individuals from the " + opt + " file '" + kFile + "'.");
	}

	private void readFamilyFile(String kFamFile, String opt) {
		BufferedReader reader = FileUtil.FileOpen(kFamFile);
		String line = null;
		famSet = NewIt.newHashSet();

		try {
			while ((line = reader.readLine()) != null) {
				String[] l = line.split(ConstValues.WHITESPACE_DELIMITER);
				famSet.add(l[0]);
			}
		} catch (IOException e) {
			Logger.handleException(e, "An exception occurred when reading the " + opt + " file '" + kFamFile + "'.");
		}
		Logger.printUserLog("Read " + famSet.size() + " family ID(s) from the " + opt + " file '" + kFamFile + "'.");
	}
}
