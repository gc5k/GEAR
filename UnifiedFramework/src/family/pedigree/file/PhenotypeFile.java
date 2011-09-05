package family.pedigree.file;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Enumeration;
import java.util.Hashtable;

import admixture.parameter.Parameter;

import family.pedigree.phenotype.FamilyUnit;
import family.pedigree.phenotype.Subject;

import util.NewIt;

/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class PhenotypeFile {

	private File phenoFile;
	private Hashtable<String, FamilyUnit> families;
	private String[] traits;
	public PhenotypeFile() {
		families = NewIt.newHashtable();
	}

	public File getPhenoFile() {
		return phenoFile;
	}

	public void parsePhenotype(File infile) throws IOException {

		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(infile));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		phenoFile = infile;
		String line;
		
		Subject sub;
		int k = 0;
		int colNum = 0;
		while ((line = reader.readLine()) != null) {
			if (line.length() == 0) {
				// skip blank lines
				continue;
			}
			if (line.startsWith("#")) {
				// skip comments
				continue;
			}
			
			if(k == 0) {
				String[] t = line.split("\\s+");
				colNum = t.length - 2;
				traits = new String[colNum];
				System.arraycopy(t, 2, traits, 0, traits.length);
				k++;
				continue;
			}
			
			String[] tokenizer = line.split("\\s+");
			int numTokens = tokenizer.length;

			sub = new Subject(colNum);
			if (numTokens < 2) {
				throw new IOException("Incorrect number of fields in phefile. line " + (k + 1));
			}

			if (colNum != numTokens - 2) {
				// this line has a different number of columns
				// should send some sort of error message
				throw new IOException("Column number mismatch in phefile. line " + (k + 1));
			}
			sub.setFamilyID(tokenizer[0]);
			sub.setSubjectID(tokenizer[1]);
			for (int j = 2; j < tokenizer.length; j++) {
				try {
					String trait = tokenizer[j];
					sub.AddTrait(trait);
				} catch (NumberFormatException nfe) {
					throw new IOException("Phenotype file input error: invalid Phenotype on line " + (k + 2));
				}
			}

			// check if the family exists already in the Hashtable
			FamilyUnit fam = this.families.get(sub.getFamilyID());
			if (fam == null) {
				// it doesnt exist, so create a new Family object
				fam = new FamilyUnit(sub.getFamilyID());
			}

			if (fam.getSubjects().containsKey(sub.getSubjectID())) {
				throw new IOException("Individual " + sub.getSubjectID() + " in family " + sub.getFamilyID() + " appears more than once.");
			}

			fam.addSubject(sub);
			this.families.put(sub.getFamilyID(), fam);
		}
		Parameter.findCovar_Number(traits);
		
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

	public int getNumTraits() {
		return traits.length;
	}

	public boolean containFamily(String fid) {
		if(families.containsKey(fid)) {
			return true;
		} else {
			return false;
		}
	}

}