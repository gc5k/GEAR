package family.qc.rowqc;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.logging.Level;

import parameter.Parameter;

import family.pedigree.Hukou;
import family.pedigree.PersonIndex;
import family.pedigree.file.MapFile;
import family.pedigree.file.PedigreeFile;
import family.pedigree.genotype.BFamilyStruct;
import family.pedigree.genotype.BPerson;
import gear.util.FileProcessor;
import gear.util.Logger;
import gear.util.NewIt;
import gear.util.TextHelper;


/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */

public class SampleFilter {

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

	public SampleFilter(PedigreeFile ped, MapFile map) {
		PedData = ped;
		MapData = map;
		PersonTable = NewIt.newArrayList();

		QC();

	}

	public void QC() {
		qualification();
		if (PersonTable.size() == 0) {
			Logger.printUserError("No individual is left for analysis.");
			System.exit(1);
		}
	}

	private void qualification() {
		ArrayList<Hukou> hukoubook = PedData.getHukouBook();
		HukouBook = NewIt.newArrayList();
		Hashtable<String, BFamilyStruct> Fam = PedData.getFamilyStruct();
		num_qualified = new int[Fam.size()][2];
		filter = new boolean[Fam.size()][];

		for (Iterator<Hukou> e = hukoubook.iterator(); e.hasNext();) {
			Hukou hukou = e.next();
			BFamilyStruct fs = Fam.get(hukou.getFamilyID());
			BPerson per = fs.getPerson(hukou.getIndividualID());
			boolean hf = hardFilter(per);
			hukou.setAvailable(hf);
			if (!hf) {
				continue;
			}
			boolean isFounder = fs.hasAncestor(per);
			HukouBook.add(hukou);
			PersonTable.add(new PersonIndex(per.getFamilyID(), per.getPersonID(), per, false, isFounder));
		}

		num_qualified = new int[Fam.size()][2];
		filter = new boolean[Fam.size()][];
		int c = 0;
		for (String fi : PedData.getFamListSorted()) {
			BFamilyStruct fs = Fam.get(fi);
			String[] pi = fs.getPersonListSorted();
			filter[c] = new boolean[pi.length];

			// filter_family

			int cc = 0;
			for (int i = 0; i < pi.length; i++) {
				BPerson per = fs.getPerson(pi[i]);
				boolean hf = hardFilter(per);
				if (!hf) {
					filter[c][cc++] = hf;
					continue;
				}

				filter[c][cc++] = hf;
				boolean isFounder = fs.hasAncestor(per);
				if (isFounder) {
					num_qualified[c][1]++;
				} else {
					num_qualified[c][0]++;
				}
			}
			c++;
		}

		for (int i = 0; i < num_qualified.length; i++) {
			qualified_Unrelated += num_qualified[i][0];
			qualified_Sib += num_qualified[i][1];
		}
	}

	protected boolean hardFilter(BPerson p) {
		boolean flag = true;

		if (Parameter.INSTANCE.keepFlag) {
			Logger.printUserLog("Reading kept individuals from '" + Parameter.INSTANCE.keepFile + "'.");
			readKeepFile();
			flag = false;
			String fi = p.getFamilyID();
			String pi = p.getPersonID();
			for (int i = 0; i < indKeep[0].length; i++) {
				if (fi.compareTo(indKeep[0][i]) == 0
						&& pi.compareTo(indKeep[1][i]) == 0) {
					flag = true;
					break;
				}
			}
		} else if (Parameter.INSTANCE.removeFlag) {
			Logger.printUserLog("Reading removed individuals from '" + Parameter.INSTANCE.removeFile + "'.");
			readRemoveFile();
			String fi = p.getFamilyID();
			String pi = p.getPersonID();
			for (int i = 0; i < indExclude[0].length; i++) {
				if (fi.compareTo(indExclude[0][i]) == 0
						&& pi.compareTo(indExclude[1][i]) == 0) {
					flag = false;
					break;
				}
			}
		}

		if (flag && Parameter.INSTANCE.keep_maleFlag) {
			return flag = p.getGender() == 1 ? true : false;
		}
		if (flag && Parameter.INSTANCE.keep_femaleFlag) {
			return flag = p.getGender() == 2 ? true : false;
		}
		if (flag && Parameter.INSTANCE.ex_nosexFlag) {
			return flag = (p.getGender() == 1 || p.getGender() == 2) ? true
					: false;
		}
		return flag;
	}

	public int getNumberMarker() {
		return MapData.getMarkerNumber();
	}

	public MapFile getMapFile() {
		return MapData;
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
		BufferedReader reader = FileProcessor.FileOpen(Parameter.INSTANCE.keepFile);
		String line = null;
		ArrayList<String> famList = NewIt.newArrayList();
		ArrayList<String> indList = NewIt.newArrayList();

		try {
			while ((line = reader.readLine()) != null) {
				String[] l = line.split(TextHelper.WHITESPACE_DELIMITER);
				if(l.length < 2) continue;
				famList.add(l[0]);
				indList.add(l[1]);
			}
		} catch (IOException e) {
			Logger.printUserError("An exception occurred when reading the kept-individual file '" + Parameter.INSTANCE.keepFile + "'.");
			Logger.printUserError("Exception Message: " + e.getMessage());
			Logger.getDevLogger().log(Level.SEVERE, "Reading kept-individual file", e);
			System.exit(1);
		}
		indKeep = new String[2][];
		indKeep[0] = (String[]) famList.toArray(new String[0]);
		indKeep[1] = (String[]) indList.toArray(new String[0]);
	}

	private void readRemoveFile() {
		BufferedReader reader = FileProcessor.FileOpen(Parameter.INSTANCE.removeFile);
		ArrayList<String> famList = NewIt.newArrayList();
		ArrayList<String> indList = NewIt.newArrayList();
		String line = null;
		try {
			while ((line = reader.readLine()) != null) {
				String[] l = line.split(TextHelper.WHITESPACE_DELIMITER);
				if(l.length < 2) continue;
				famList.add(l[0]);
				indList.add(l[1]);
			}
		} catch (IOException e) {
			Logger.printUserError("An exception occurred when reading the removed-individual file '" + Parameter.INSTANCE.removeFile + "'.");
			Logger.printUserError("Exception Message: " + e.getMessage());
			Logger.getDevLogger().log(Level.SEVERE, "Reading removed-individual file", e);
			System.exit(1);
		}
		indExclude = new String[2][];
		indExclude[0] = (String[]) famList.toArray(new String[0]);
		indExclude[1] = (String[]) indList.toArray(new String[0]);

	}
}
