package family.pedigree.file;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.Hashtable;

import admixture.parameter.Parameter;

import util.NewIt;
import family.mdr.MDRConstant;
import family.pedigree.genotype.BFamilyStruct;
import family.pedigree.genotype.BPerson;

/**
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class PedigreeFile {

	private ArrayList<String[]> AlleleSet;
	private Hashtable<String, Boolean> famInformative;
	private Hashtable<String, BFamilyStruct> familystructure;

	private boolean bogusParents = false;
	// private ArrayList<SNP> markerInfor;
	private int num_marker;
	private ArrayList<String> pedigrees;
	private String titleLine = null;
	private File pedfile;
	private boolean header = true;
	static int MissingAllele = 0;

	// public static String MissingGenotype="00";
	public PedigreeFile() {
		AlleleSet = NewIt.newArrayList();
		this.famInformative = NewIt.newHashtable();
		this.familystructure = NewIt.newHashtable();
	}

	public String[] getFamListSorted() {
		Enumeration<String> famstrList = this.familystructure.keys();
		String[] FID = new String[familystructure.size()];
		int ind = 0;
		while (famstrList.hasMoreElements()) {
			FID[ind++] = (String) famstrList.nextElement();
		}
		Arrays.sort(FID);
		return FID;
	}

	public int getNumInformativeFamilies() {
		int c = 0;
		Enumeration<String> famstrList = familystructure.keys();
		while (famstrList.hasMoreElements()) {
			c += famInformative.get(famstrList.nextElement()).booleanValue() ? 1 : 0;
		}
		return c;
	}

	public String[] getInformativeFamListSorted() {
		Enumeration<String> famstrList = this.familystructure.keys();
		String[] FID = new String[getNumInformativeFamilies()];
		int ind = 0;
		while (famstrList.hasMoreElements()) {
			String F = (String) famstrList.nextElement();
			if (famInformative.get(F).booleanValue()) {
				FID[ind++] = F;
			}
		}
		Arrays.sort(FID);
		return FID;
	}

	public String[] getUninformativeFamListSorted() {
		Enumeration<String> famstrList = this.familystructure.keys();
		String[] FID = new String[familystructure.size() - getNumInformativeFamilies()];
		int ind = 0;
		while (famstrList.hasMoreElements()) {
			String F = (String) famstrList.nextElement();
			if (!famInformative.get(F).booleanValue()) {
				FID[ind++] = F;
			}
		}
		Arrays.sort(FID);
		return FID;
	}

	public Hashtable<String, Boolean> getFamInformative() {
		return this.famInformative;
	}

	public BFamilyStruct getFamilyStruct(String familystrID) {
		return this.familystructure.get(familystrID);
	}

	/**
	 * this method iterates through each family in Hashtable families and adds
	 * up the number of individuals in total across all families
	 * 
	 * @return the total number of individuals in all the family objects in the
	 *         families hashtable
	 */
	public int getNumIndividuals() {
		Enumeration<BFamilyStruct> famEnum = familystructure.elements();
		int total = 0;
		while (famEnum.hasMoreElements()) {
			BFamilyStruct fam = famEnum.nextElement();
			total += fam.getNumPersons();
		}
		return total;
	}

	public Hashtable<String, BFamilyStruct> getFamilyStruct() {
		return familystructure;
	}

	public File GetPedFile() {
		return pedfile;
	}

	public void Initial(File infile) throws MDRPedFileException, IOException {
		pedigrees = NewIt.newArrayList();
		BufferedReader reader = new BufferedReader(new FileReader(infile));
		pedfile = infile;
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
			pedigrees.add(line);
		}
		if (header) {
			titleLine = (String) pedigrees.get(0);
			pedigrees.remove(0);

			int numLines = pedigrees.size();

			if (numLines < 2) {
				throw new MDRPedFileException("Pedgree data format error: empty pedigree file");
			}
			String[] tokenizer = titleLine.split("\\s+");
			int numTokens = tokenizer.length;

			if (numTokens < 7) {
				throw new MDRPedFileException("Pedgree format error: the title line is incorrect");
			}
			num_marker = numTokens - 6;

			// System.out.println(markerInfor);
		} else {
			int numLines = pedigrees.size();
			if (numLines < 1) {
				throw new MDRPedFileException("Pedgree data format error: empty pedigree file");
			}
			String[] tokenizer = pedigrees.get(0).split("\\s+");
			num_marker = (tokenizer.length - 6) / 2;
		}
		AlleleSet.ensureCapacity(num_marker);
		for (int i = 0; i < num_marker; i++) {
			String[] a = new String[2];
			a[0] = a[1] = Parameter.missing_allele;
			AlleleSet.add(a);
		}
	}

	/**
	 * Taking in a pedigree file in the form of a vector of strings and parses
	 * it. The data parsed is stored in families in the member hashtable
	 * families. Note that the "Linkage" here is the relationship between
	 * relatives in a pedigree, but is not that term of genetics.
	 */
	public void parseLinkage() throws MDRPedFileException {
		int colNum = num_marker * 2 + 6;

		int numMarkers = 0;
		int numLines = pedigrees.size();
		if (numLines == 0) {
			throw new MDRPedFileException("Data format error: empty file");
		}
		BPerson per;

		for (int k = 0; k < numLines; k++) {
			String[] tokenizer = pedigrees.get(k).split(MDRConstant.delim);

			int numTokens = tokenizer.length;

			if (numTokens % 2 == 0) {
				numMarkers = (numTokens - 6) / 2;
				if (numMarkers != num_marker) {
					new MDRPedFileException("Mismatched Colunm in pedfile. line " + (k + 2));
				}
			} else {
				new MDRPedFileException("Mismatched Column in pedfile. line " + (k + 2));
			}

			if (colNum != numTokens) {
				// this line has a different number of columns
				// should send some sort of error message
				throw new MDRPedFileException("Column number mismatch in pedfile. line " + (k + 2));
			}

			per = new BPerson(numMarkers);
			// pseudoper = new PseudoPerson();
			if (numTokens < 6) {
				throw new MDRPedFileException("Incorrect number of fields on line " + (k + 2));
			}

			if (tokenizer.length > 6) {

				per.setFamilyID(tokenizer[0]);
				per.setPersonID(tokenizer[1]);
				per.setDadID(tokenizer[2]);
				per.setMomID(tokenizer[3]);

				try {
					int Gender = Integer.parseInt(tokenizer[4]);
					int Status = Integer.parseInt(tokenizer[5]);

					per.setGender(Gender);
					per.setAffectedStatus(Status);

				} catch (NumberFormatException nfe) {
					throw new MDRPedFileException("Pedfile error: invalid gender or affected status on line " + (k + 2));
				}

				for (int j = 0; j < (tokenizer.length - 6) / 2; j++) {
					try {
						String[] allele = {tokenizer[6 + j * 2], tokenizer[6 + j * 2 + 1]};
						boolean flag = (allele[0].compareTo(Parameter.missing_allele) != 0) && (allele[1].compareTo(Parameter.missing_allele) != 0);
						if (flag) {
							int[] code = recode(j, allele);
							per.addMarker(flag, code[0], code[1], j);
						} else {
							per.addMarker(flag, 0, 0, j);
						}

					} catch (NumberFormatException nfe) {
						throw new MDRPedFileException("Pedigree file input error: invalid genotype on line " + (k + 2));
					}
				}
				// check if the family exists already in the Hashtable
				BFamilyStruct famstr = familystructure.get(per.getFamilyID());
				if (famstr == null) {
					// it doesn't exist, so create a new FamilyStruct object
					famstr = new BFamilyStruct(per.getFamilyID());
					familystructure.put(per.getFamilyID(), famstr);
				}

				if (famstr.getPersons().containsKey(per.getPersonID())) {
					throw new MDRPedFileException("Person " + per.getPersonID() + " in family " + per.getFamilyID() + " appears more than once.");
				}
				famstr.addPerson(per);
				if (Parameter.mode.compareTo("pi") == 0) {
					famstr.addPseudoPerson(per);
				}
			}
		}

		pedigrees.clear();
		pedigrees = null;
	}

	private int[] recode(int idx, String[] allele) {
		int[] code = { -1, -1 };
		String[] ref = AlleleSet.get(idx);
		if (ref[1] != Parameter.missing_allele) { // two detected alleles
			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < 2; j++) {
					if (allele[i].compareTo(ref[j]) == 0) {
						code[i] = j;
						break;
					}
				}
			}
			if (code[0] == -1 || code[1] == -1) {
				System.err.println("more than 2 alleles in the marker column " + (idx + 1));
			}
		} else {
			// less than two detected alleles
			
			if (allele[0].compareTo(allele[1]) == 0) {
				// 1 when both alleles are same
				if (ref[0].compareTo(Parameter.missing_allele)==0) {
					// zero detected alleles
					ref[0] = new String(allele[0]);
					code[0] = code[1] = 0;
				} else {
					// one detected alleles
					if (allele[0].compareTo(ref[0]) == 0) {
						code[0] = code[1] = 0;
					} else {
						code[0] = code[1] = 1;
						ref[1] = new String(allele[1]);
					}
				}
			} else {
				// 2 when both alleles are different
				if (ref[0].compareTo(Parameter.missing_allele) == 0) {
					// zero detected alleles
					ref[0] = new String(allele[0]); ref[1] = new String(allele[1]);
					code[0] = 0; code[1] = 1;
				} else {
					// one detected alleles
					if (ref[0].compareTo(allele[0]) == 0) {
						ref[1] = new String(allele[1]);
						code[0] = 0; code[1] = 1;
					} else if (ref[0].compareTo(allele[1]) == 0) {
						ref[1] = new String(allele[0]);
						code[0] = 1; code[1] = 0;
					} else {
						System.out.println("more than 3 alleles in marker column " + (idx + 1));
					}
				}
			}
		}
		return code;
	}

	public ArrayList<String[]> getPolymorphism() {
		return AlleleSet;
	}

	public boolean hasBogusParents() {
		return bogusParents;
	}

	public int getNumMarker() {
		return num_marker;
	}

	public void setHeader(boolean flag) {
		header = flag;
	}
}