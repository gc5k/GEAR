package family.pedigree.file;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Hashtable;

import family.pedigree.phenotype.FamilyUnit;
import family.pedigree.phenotype.Subject;

import util.NewIt;

/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class GMDRPhenoFile {

	private String titleLine;
	private ArrayList<String> traitInfor;
	private ArrayList<String> traitStrings;
	private File phenoFile;
	private ArrayList<Subject> allSubjects;
	private Hashtable<String, FamilyUnit> families;

	public GMDRPhenoFile() {
		families = NewIt.newHashtable();
	}

	public void Initial(File infile) throws GMDRPhenoFileException, IOException {
		traitStrings = NewIt.newArrayList();
		BufferedReader reader = new BufferedReader(new FileReader(infile));
		phenoFile = infile;
		String line;
		while ((line = reader.readLine()) != null) {
			if (line.length() == 0) {
				// skip blank lines
				continue;
			}
			if (line.startsWith("#")) {
				// skip comments
				continue;
			}
			traitStrings.add(line);
		}
		titleLine = (String) traitStrings.get(0);
		traitStrings.remove(0);
		int numLines = traitStrings.size();

		if (numLines < 2) {
			throw new GMDRPhenoFileException("Phenotype data format error: empty file");
		}
		String[] tokenizer = titleLine.split("\\s+");

		int numTokens = tokenizer.length;

		if (numTokens < 2) {
			throw new GMDRPhenoFileException("Phenotype data format error: the title line is incorrect");
		}

		// reading the title line:get the marker number
		traitInfor = NewIt.newArrayList();
		for(int i = 2; i < tokenizer.length; i++) {
			traitInfor.add(tokenizer[i]);
		}
	}

	public File getPhenoFile() {
		return phenoFile;
	}

	public void parsePhenotype() throws GMDRPhenoFileException {
		int colNum = -1;
		if (traitStrings.size() < 2) {
			throw new GMDRPhenoFileException("Pheno Data format error: empty file");
		}
		Subject sub;
		this.allSubjects = NewIt.newArrayList();

		for (int k = 0; k < traitStrings.size(); k++) {
			String[] tokenizer = traitStrings.get(k).split("\\s+");
			int numTokens = tokenizer.length;

			sub = new Subject(traitInfor.size());
			if (colNum < 1) {
				// only check column number count for the first nonblank line
				colNum = numTokens;
			}
			if (numTokens < 2) {
				throw new GMDRPhenoFileException("Incorrect number of fields in phefile. line " + (k + 1));
			}

			if (colNum != numTokens || colNum != traitInfor.size() + 2) {
				// this line has a different number of columns
				// should send some sort of error message
				throw new GMDRPhenoFileException("Column number mismatch in phefile. line " + (k + 2));
			}
			sub.setFamilyID(tokenizer[0]);
			sub.setSubjectID(tokenizer[1]);
			for (int j = 2; j < tokenizer.length; j++) {
				try {
					String trait = tokenizer[j];
					sub.AddTrait(trait);
				} catch (NumberFormatException nfe) {
					throw new GMDRPhenoFileException("Phenotype file input error: invalid Phenotype on line " + (k + 2));
				}
			}

			// check if the family exists already in the Hashtable
			FamilyUnit fam = this.families.get(sub.getFamilyID());
			if (fam == null) {
				// it doesnt exist, so create a new Family object
				fam = new FamilyUnit(sub.getFamilyID());
			}

			if (fam.getSubjects().containsKey(sub.getSubjectID())) {
				throw new GMDRPhenoFileException("Individual " + sub.getSubjectID() + " in family " + sub.getFamilyID() + " appears more than once.");
			}

			fam.addSubject(sub);
			this.families.put(sub.getFamilyID(), fam);
			this.allSubjects.add(sub);
		}
		traitStrings = null;
	}

	public FamilyUnit getFamilyUnit(String familyUnitID) {
		return this.families.get(familyUnitID);
	}

	public Enumeration<String> getFamUnitList() {
		return this.families.keys();
	}

	public int getNumFamilyUnits() {
		return families.size();
	}

	public int getNumSubjects() {
		return allSubjects.size();
	}

	public int getNumTraits() {
		return traitInfor.size();
	}

	public ArrayList<String> getTraitName() {
		return traitInfor;
	}

	public boolean containFamily(String fid) {
		if(families.containsKey(fid)) {
			return true;
		} else {
			return false;
		}
	}

	public String getTraitAtI(int index) {
		if (index < 0 || index > traitInfor.size()) {
			System.err.println("Could not find the phenotype at index " + index);
			System.exit(0);
		}
		return traitInfor.get(index);
	}
}