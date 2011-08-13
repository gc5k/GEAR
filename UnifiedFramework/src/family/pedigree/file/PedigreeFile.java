package family.pedigree.file;

/*
 * $Id: PedFile.java,v 3.21 2006/05/12 18:01:28 jmaller Exp $
 * WHITEHEAD INSTITUTE
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2002 by the
 * Whitehead Institute for Biomedical Research.  All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support
 * whatsoever.  The Whitehead Institute can not be responsible for its
 * use, misuse, or functionality.
 */
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.Hashtable;

import admixture.parameter.Parameter;

import publicAccess.PublicData;
import util.NewIt;
import family.pedigree.genotype.FamilyStruct;
import family.pedigree.genotype.Person;

/**
 * Handles input and storage of Pedigree files this class is not thread safe
 * (untested). modified from original Pedfile and checkdata classes by Hui Gong
 * 
 * @author Julian Maller
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class PedigreeFile {

	private ArrayList<HashSet<String>> AlleleSet;
	private Hashtable<String, Boolean> famInformative;
	private Hashtable<String, FamilyStruct> familystructure;

	// stores the individuals found by parse() in allIndividuals. this is useful
	// for outputting Pedigree information to
	// a file of another type.

	// stores the individuals chosen by pedparser
	// bogusParents is true if someone in the file referenced a parent not in
	// the file
	private boolean bogusParents = false;
//	private ArrayList<SNP> markerInfor;
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

	public FamilyStruct getFamilyStruct(String familystrID) {
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
		Enumeration<FamilyStruct> famEnum = familystructure.elements();
		int total = 0;
		while (famEnum.hasMoreElements()) {
			FamilyStruct fam = famEnum.nextElement();
			total += fam.getNumPersons();
		}
		return total;
	}

	public Hashtable<String, FamilyStruct> getFamilyStruct() {
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
			num_marker = (tokenizer.length - 6)/2;
		}
		AlleleSet.ensureCapacity(num_marker);
		for(int i = 0; i < num_marker; i++) {
			AlleleSet.add(new HashSet<String>());
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
		boolean genoError = false;
		int numLines = pedigrees.size();
		if (numLines == 0) {
			throw new MDRPedFileException("Data format error: empty file");
		}
		Person per;

		for (int k = 0; k < numLines; k++) {
			String[] tokenizer = pedigrees.get(k).split(PublicData.delim);

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

			per = new Person(numMarkers);
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

				byte genotype1;
				byte genotype2;
				for (int j = 0; j < (tokenizer.length - 6) / 2; j++) {
					try {
						String alleleA = tokenizer[6 + j * 2];
						String alleleB = tokenizer[6 + j * 2 + 1];
						int[] checker1, checker2;
						checker1 = checkGenotype(alleleA);
						checker2 = checkGenotype(alleleB);
						if (checker1[1] != checker2[1]) {
							genoError = !genoError;
						}
						Polymorphism(j, alleleA, alleleB);
						if (genoError) {
							throw new MDRPedFileException("File input error on line " + (k + 2) + ", marker " + (per.getNumMarkers() + 2)
									+ ".\nFor any marker, an individual's genotype must be only letters or only numbers.");
						}
						if (checker1[0] == PublicData.MissingAllele || checker2[0] == PublicData.MissingAllele) {
							checker1[0] = checker2[0] = PublicData.MissingAllele;
						}
						genotype1 = (byte) checker1[0];
						genotype2 = (byte) checker2[0];
						per.addMarker(genotype1, genotype2);
					} catch (NumberFormatException nfe) {
						throw new MDRPedFileException("Pedigree file input error: invalid genotype on line " + (k + 2));
					}
				}
				// check if the family exists already in the Hashtable
				FamilyStruct famstr = familystructure.get(per.getFamilyID());
				if (famstr == null) {
					// it doesn't exist, so create a new FamilyStruct object
					famstr = new FamilyStruct(per.getFamilyID());
					familystructure.put(per.getFamilyID(), famstr);
				}

				if (famstr.getPersons().containsKey(per.getPersonID())) {
					throw new MDRPedFileException("Person " + per.getPersonID() + " in family " + per.getFamilyID() + " appears more than once.");
				}

				famstr.addPerson(per);
				// famstr.addPseudoPerson(pseudoper);
			}
		}

		pedigrees.clear();
		pedigrees = null;
	}

	public int[] checkGenotype(String allele) throws MDRPedFileException {
		// This method cleans up the genotype checking process for hap map and
		// ped files & allows for both numerical and
		// alphabetical input.
		int[] genotype = new int[2];

		if (allele.compareTo(Parameter.missing_allele) == 0) {
			genotype[0] = 0;
		} else if (allele.equalsIgnoreCase("A")) {
			genotype[0] = 1;
		} else if (allele.equalsIgnoreCase("C")) {
			genotype[0] = 2;
		} else if (allele.equalsIgnoreCase("G")) {
			genotype[0] = 3;
		} else if (allele.equalsIgnoreCase("T")) {
			genotype[0] = 4;
		} else {
			genotype[0] = Integer.parseInt(allele.trim());
			genotype[1] = 1;
		}

		return genotype;
	}

	private void Polymorphism(int i, String A, String B) {
		if(A.compareTo(Parameter.missing_allele) != 0) {
			AlleleSet.get(i).add(A);
		} 
		if(B.compareTo(Parameter.missing_allele) != 0){
			AlleleSet.get(i).add(B);
		}
	}

	public ArrayList<HashSet<String>> getPolymorphism() {
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