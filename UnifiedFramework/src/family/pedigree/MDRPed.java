package family.pedigree;

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
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.TreeMap;

import publicAccess.PublicData;
import util.NewIt;
import family.RabinowitzLairdAlgorithm.AbstractGenoDistribution;
import family.imputation.*;
import family.pedigree.genotype.FamilyStruct;
import family.pedigree.genotype.FamilyStructException;
import family.pedigree.genotype.Person;
import family.pedigree.genotype.PseudoPerson;
import edu.mit.wi.haploview.Chromosome;
import edu.mit.wi.haploview.SNP;
import edu.mit.wi.pedfile.MarkerResult;


/**
 * Handles input and storage of Pedigree files this class is not thread safe
 * (untested). modified from original Pedfile and checkdata classes by Hui Gong
 * 
 * @author Julian Maller
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class MDRPed {

	private Hashtable<String, Boolean> famInformative;
	private Hashtable<String, FamilyStruct> familystructure;
	private ArrayList<Person> axedPeople = NewIt.newArrayList();

	// stores the individuals found by parse() in allIndividuals. this is useful
	// for outputting Pedigree information to
	// a file of another type.
	private ArrayList<Person> allIndividuals;

	// stores the individuals chosen by pedparser
	private ArrayList<MarkerResult> results = null;
	private String[][] hminfo;
	// bogusParents is true if someone in the file referenced a parent not in
	// the file
	private boolean bogusParents = false;
	private int[] markerRatings;
	private int[] dups;
	private ArrayList<String> markerInfor;
	private ArrayList<String> pedigrees;
	private ArrayList<MendErrorTrace> menderrortrace;
	private String titleLine;
	private File pedfile;
	static int MissingAllele = 0;

	// public static String MissingGenotype="00";
	public MDRPed() {

		// hardcoded hapmap info
		this.famInformative = NewIt.newHashtable();
		this.familystructure = NewIt.newHashtable();
		this.menderrortrace = NewIt.newArrayList();
	}

	public void GenotypeImputation() throws MDRPedFileException {
		Enumeration<String> fsList = this.familystructure.keys();
		while (fsList.hasMoreElements()) {
			FamilyStruct fs = (FamilyStruct) getFamilyStruct((String) fsList.nextElement());
			for (int i = 0; i < getNumMarkers(); i++) {
				System.out.println(fs.getFamilyStructName() + " " + i);
				Imputation(fs, i);
			}
		}
		GenotypeSummary(false);
	}

	protected void Imputation(FamilyStruct fs, int genoIdx) throws MDRPedFileException {
		GenoSet gSet = fs.getObservedGenoSet(genoIdx);
		AbstractImputation ai;
		if ((gSet.getNumUntypedChildren() > 0) && (gSet.getNumUntypedChildren() < gSet.getNumChildren())) {
			if (gSet.getNumTypedParents() == 2) {
				ai = new GenotypedParents(gSet.getchildrenGenoMap(), gSet.getparentsGenoMap());
			} else if (gSet.getNumParents() == 1) {
				String pg = gSet.getfullchildrenGenoMap().firstKey();
				if (!AbstractGenoDistribution.isHeterozygous(pg)) {
					ai = new OneHomozygousParent(gSet.getchildrenGenoMap(), gSet.getparentsGenoMap());
				} else {
					ai = new OneHeterozygousParent(gSet.getchildrenGenoMap(), gSet.getparentsGenoMap());
				}
			} else {
				ai = new UngenotypedParents(gSet.getchildrenGenoMap());
			}

			Enumeration<String> perList = fs.getPersonList();
			while (perList.hasMoreElements()) {
				Person per;
				per = fs.getPerson(perList.nextElement());
				String genotype = per.getGenotype(genoIdx);
				if (fs.hasAncestor(per.getPersonID())) {
					if (genotype.compareTo(PublicData.MissingGenotype) == 0) {
						String ImputedGenotype = ai.RandomAssign();
						per.setGenotype(genoIdx, ImputedGenotype);
						System.out.println(fs.getFamilyStructName() + " " + per.getPersonID() + " " + genoIdx);
					}
				}
			}
		}
	}

	protected void GenotypeSummary(boolean IsObservedGenoSet) throws MDRPedFileException {
		Enumeration<String> fsList = this.familystructure.keys();
		ArrayList<String> error_report = NewIt.newArrayList();
		while (fsList.hasMoreElements()) {
			FamilyStruct fs = (FamilyStruct) getFamilyStruct((String) fsList.nextElement());
			Boolean b = new Boolean(true);
			famInformative.put(fs.getFamilyStructName(), b);
			TreeMap<String, Integer> Ps;
			TreeMap<String, Integer> Ks;
			GenoSet gSet;
			for (int i = 0; i < getNumMarkers(); i++) {
				Ps = NewIt.newTreeMap();
				Ks = NewIt.newTreeMap();
				Enumeration<String> perList = fs.getPersonList();
				while (perList.hasMoreElements()) {
					Person per = fs.getPerson(perList.nextElement());
					String genotype = per.getGenotype(i);
					if (fs.hasAncestor(per.getPersonID())) {
						if (Ks.containsKey(genotype)) {
							Integer c = ((Integer) Ks.get(genotype));
							c++;
							Ks.put(genotype, c);
						} else {
							Integer c = new Integer(1);
							Ks.put(new String(genotype), c);
						}
					} else {
						if (Ps.containsKey(genotype)) {
							Integer c = ((Integer) Ps.get(genotype));
							c++;
							Ps.put(genotype, c);
						} else {
							Integer c = new Integer(1);
							Ps.put(new String(genotype), c);
						}
					}
				}
				gSet = new GenoSet(Ps, Ks, i);
				if (IsObservedGenoSet) {
					fs.addObservedGenoSet(gSet);
				} else {
					if (gSet.getchildrenGenoMap().size() == 0) {
						error_report.add(new String("Family " + fs.getFamilyStructName() + " Locus " + markerInfor.get(i)));
						b = new Boolean(false);
						famInformative.put(fs.getFamilyStructName(), b);
					}
					fs.addImputedGenoSet(gSet);
				}
			}
		}
	}

	public Enumeration<String> getFamStrList() {
		return this.familystructure.keys();
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
	/**
	 * finds the first individual in the first family and returns the number of
	 * markers for that individual
	 * 
	 * @return the number of markers
	 */
	public int getNumMarkers() {
		Enumeration<FamilyStruct> famList = familystructure.elements();
		int numMarkers = 0;
		while (famList.hasMoreElements()) {
			FamilyStruct fam = (FamilyStruct) famList.nextElement();
			Enumeration indList = fam.getPersonList();
			Person per = null;
			while (indList.hasMoreElements()) {
				per = fam.getPerson((String) indList.nextElement());
				numMarkers = per.getNumMarkers();
				if (numMarkers > 0) {
					return numMarkers;
				}
			}
		}
		return 0;
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
		titleLine = (String) pedigrees.get(0);
		pedigrees.remove(0);
		int numLines = pedigrees.size();

		if (numLines < 2) {
			throw new MDRPedFileException("Pedgree data format error: empty file");
		}
		StringTokenizer tokenizer = new StringTokenizer(titleLine, "\n\t\" \"");
		int numTokens = tokenizer.countTokens();

		if (numTokens < 7) {
			throw new MDRPedFileException("Pedgree format error: the title line is incorrect");
		}

		// reading the title line:get the marker number
		int c = 0;
		ArrayList<String> temp = NewIt.newArrayList();
		while (tokenizer.hasMoreTokens()) {
			if (c++ < 6) {
				tokenizer.nextToken();
			} else {
				String marker = (String) tokenizer.nextToken();
				temp.add(marker);
			}
		}
		markerInfor = temp;
		// System.out.println(markerInfor);
		
	}

	/**
	 * Taking in a pedigree file in the form of a vector of strings and parses
	 * it. The data parsed is stored in families in the member hashtable
	 * families. Note that the "Linkage" here is the relationship between
	 * relatives in a pedigree, but is not that term of genetics.
	 */
	public void parseLinkage() throws MDRPedFileException {
		int colNum = markerInfor.size() * 2 + 6;

		int numMarkers = 0;
		boolean genoError = false;
		int numLines = pedigrees.size();
		if (numLines == 0) {
			throw new MDRPedFileException("Data format error: empty file");
		}
		Person per;
		PseudoPerson pseudoper;
		this.allIndividuals = NewIt.newArrayList();

		for (int k = 0; k < numLines; k++) {
			String[] tokenizer = pedigrees.get(k).split(PublicData.delim);

			int numTokens = tokenizer.length;

			if (numTokens % 2 == 0) {
				numMarkers = (numTokens - 6) / 2;
				if (numMarkers != markerInfor.size()) {
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
			pseudoper = new PseudoPerson();
			if (numTokens < 6) {
				throw new MDRPedFileException("Incorrect number of fields on line " + (k + 2));
			}

			if (tokenizer.length > 6) {

				per.setFamilyID(tokenizer[0]);
				per.setPersonID(tokenizer[1]);
				per.setDadID(tokenizer[2]);
				per.setMomID(tokenizer[3]);

				pseudoper.setFamilyID(tokenizer[0]);
				pseudoper.setPseudoPersonID(tokenizer[1]);
				pseudoper.setDadID(tokenizer[2]);
				pseudoper.setMomID(tokenizer[3]);

				try {
					int Gender = Integer.parseInt(tokenizer[4]);
					int Status = Integer.parseInt(tokenizer[5]);

					per.setGender(Gender);
					per.setAffectedStatus(Status);

					pseudoper.setGender(Gender);
					pseudoper.setAffectedStatus(Status);
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
					famstr.setNumMarkers(numMarkers);
					familystructure.put(per.getFamilyID(), famstr);
				}

				if (famstr.getPersons().containsKey(per.getPersonID()) || famstr.getPseudoPersons().containsKey(pseudoper.getPseudoPersonID())) {
					throw new MDRPedFileException("Person " + per.getPersonID() + " in family " + per.getFamilyID() + " appears more than once.");
				}

				famstr.addPerson(per);
				famstr.addPseudoPerson(pseudoper);
			}
		}

		// now we check if anyone has a reference to a parent who isn't in the
		// file, and if so, we remove the reference
		for (int i = 0; i < allIndividuals.size(); i++) {
			Person currentInd = (Person) allIndividuals.get(i);
			Hashtable curFam = familystructure.get(currentInd.getFamilyID()).getPersons();
			if (!currentInd.getDadID().equals("0") && !(curFam.containsKey(currentInd.getDadID()))) {
				currentInd.setDadID("0");
				bogusParents = true;
			}
			if (!currentInd.getMomID().equals("0") && !(curFam.containsKey(currentInd.getMomID()))) {
				currentInd.setMomID("0");
				bogusParents = true;
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

		if (allele.equalsIgnoreCase("N")) {
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

	public String[][] getHMInfo() {
		return hminfo;
	}

	public ArrayList getResults() {
		return results;
	}

	public void setResults(ArrayList res) {
		results = res;
	}

	public ArrayList getAxedPeople() {
		return axedPeople;
	}

	public boolean isBogusParents() {
		return bogusParents;
	}

	public ArrayList<String> getMarkerInformation() {
		return markerInfor;
	}

	public ArrayList<String> getMarkerInformation(int[] subsetMarker) {
		if (subsetMarker.length == markerInfor.size() || subsetMarker == null) {
			return markerInfor;
		} else {
			ArrayList<String> mk = NewIt.newArrayList();
			for (int i = 0; i < subsetMarker.length; i++) {
				mk.add(markerInfor.get(subsetMarker[i]));
			}
			return mk;
		}
	}

	public int[] getMarkerRatings() {
		return markerRatings;
	}

	public int[] getDups() {
		return dups;
	}

	public void RabinowitzApproach(boolean forNontransmitted, int[] subsetMarker) {
		Enumeration fsList = this.familystructure.keys();
		boolean informative;
		String fid;
		while (fsList.hasMoreElements()) {
			fid = (String) fsList.nextElement();
			FamilyStruct fs = (FamilyStruct) getFamilyStruct(fid);
			if (!((Boolean) famInformative.get(fid)).booleanValue()) {
				// System.err.println("Omitted" + fs.getFamilyStructName());
				continue;
			}
			try {
				if (forNontransmitted) {
					fs.NontransmittedProc(markerInfor, subsetMarker);
				} else {
					fs.RabinowitzProc(markerInfor, subsetMarker);
				}
			} catch (FamilyStructException E) {
				System.err.println("Exception in family " + fs.getFamilyStructName() + " when in RabinowitzApproach.");
			}
		}
	}

	/**
	 * return the record of menderrortrace(); Allele2Genotype should be invoked
	 * precedingly;
	 * 
	 * @return
	 */
	public ArrayList getMendErrorTrace() {
		return menderrortrace;
	}

	public void printFamilyStruct() {
		String[] FID = new String[famInformative.size()];
		int idx = 0;
		Set keys = famInformative.keySet();
		for (Iterator i = keys.iterator(); i.hasNext();) {
			FID[idx++] = (String) i.next();
		}
		Arrays.sort(FID);

		for (int i = 0; i < FID.length; i++) {
			String fid = FID[i];
			FamilyStruct fs = getFamilyStruct(fid);
			Hashtable htPerson = fs.getPseudoPersons();
			Enumeration pid = htPerson.keys();
			String[] PID = new String[htPerson.size()];
			idx = 0;
			while (pid.hasMoreElements()) {
				PID[idx++] = (String) pid.nextElement();
			}
			Arrays.sort(PID);
			for (int j = 0; j < PID.length; j++) {
				PseudoPerson ps = (PseudoPerson) htPerson.get(PID[j]);
				System.out.print(ps.getPseudoPersonID() + "\t");
				ArrayList mk = (ArrayList) ps.getGenotype();
				for (int k = 0; k < mk.size(); k++) {
					System.out.print(mk.get(k) + "\t");
				}
				System.out.println();
			}
		}
	}
}