package family.pedigree.genotype;

/*
 * $Id: Family.java,v 3.1 2006/04/10 18:29:51 djbender Exp $
 * WHITEHEAD INSTITUTE
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2002 by the
 * Whitehead Institute for Biomedical Research.  All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support
 * whatsoever.  The Whitehead Institute can not be responsible for its
 * use, misuse, or functionality.
 */
import java.util.Arrays;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Set;
import java.util.TreeMap;

import util.NewIt;
import family.RabinowitzLairdAlgorithm.AbstractGenoDistribution;

/**
 * Storing the familyName and the members of a family from a pedigree file. This class is not thread safe (untested)
 * 
 * @author Julian Maller, extended by Guo-Bo Chen, chenguobo@gmail.com
 */
public class FamilyStruct {
	// save observed genotypes;
	private Hashtable<String, Person> persons;
	private Hashtable<String, PseudoPerson> pseudopersons;
	private String familyStructName;
	private int mendErrors;
	private int numMarkers;

	/**
	 * adds a member to this family (adds to persons ArrayList)
	 * 
	 * @param ind
	 *            Person to add to persons ArrayList
	 */
	public void addPerson(Person per) {
		this.persons.put(per.getPersonID(), per);
	}

	public void addPseudoPerson(PseudoPerson pseudoper) {
		this.pseudopersons.put(pseudoper.getPseudoPersonID(), pseudoper);
	}

	public void countAllele(TreeMap<String, Integer> Geno, Set<String> alleleSet) {

		for (String g : Geno.keySet()) {
			alleleSet.add(g.substring(0, 1));
			alleleSet.add(g.substring(1, 2));
		}
	}

	public boolean containsPerson(String id) {
		if (persons.containsKey(id)) {
			return true;
		} else {
			return false;
		}
	}

	public boolean containsPseudoPerson(String id) {
		if (pseudopersons.containsKey(id)) {
			return true;
		} else {
			return false;
		}
	}

	public FamilyStruct() {
		this.persons = NewIt.newHashtable();
		this.pseudopersons = NewIt.newHashtable();
	}

	public FamilyStruct(String familyStructName) {
		this.persons = NewIt.newHashtable();
		this.pseudopersons = NewIt.newHashtable();
		this.familyStructName = familyStructName;
	}

	/**
	 * returns the name of this family (familyName)
	 * 
	 * @return family name
	 */
	public String getFamilyStructName() {
		return familyStructName;
	}

	/**
	 * Table 1: Conditional distribution of non-transmitted genotypes when one homozygous parents's genotype and children's genotypes are available at
	 * marker locus i Table 2: Conditional contribution of non-transmitted genotypes when one heterozygous parent's genotype and children's genotype
	 * are available at marker locus i Table 3: Conditional contribution of non-transmitted genotypes when only children's genotypes but no parent's
	 * genotype are available at marker locus A
	 * 
	 * @param genotype
	 *            , transmitted genotype at locus i of indiviudal j
	 * @param genoset
	 *            , configuration of parents' and children's genotypes.
	 * @return non-transmitted genotype
	 */
	public String[] getNonTransmitted(String transmitted, AbstractGenoDistribution genodis) {
		String nontran_tran[] = genodis.getNontransmitted(transmitted);
		// System.out.println(nontran_tran[0] + "\t" + nontran_tran[1]);
		return nontran_tran;

	}

	public int getMendErrs() {
		return mendErrors;
	}

	/**
	 * returns the number of persons of this family
	 * 
	 * @return number of persons in this family
	 */
	public int getNumPersons() {
		return this.persons.size();
	}

	public int getNumSibs() {
		int i = 0;
		for (String pi:persons.keySet()) {
			if(hasAncestor(pi)) {
				i++;
			}
		}
		return i;
	}
	
	public int getNumFounders() {
		int i = 0; 
		for (String pi:persons.keySet()) {
			if(!hasAncestor(pi)) {
				i++;
			}
		}
		return i;
	}

	/**
	 * returns a list of personIDs that are persons of this family in the form of an enumeration
	 * 
	 * @return enumeration memberlist
	 */
	public Enumeration<String> getPersonList() {
		return this.persons.keys();
	}

	public String[] getPersonListSorted() {
		Enumeration<String> perstrList = getPersonList();
		String[] PID = new String[persons.size()];
		int ind = 0;
		while (perstrList.hasMoreElements()) {
			PID[ind++] = perstrList.nextElement();
		}
		Arrays.sort(PID);
		return PID;
	}

	/**
	 * get the Person with personID
	 * 
	 * @param personID
	 *            the personID of the person we want
	 * @return the person with matching personID
	 */
	public Person getPerson(String personID) {
		return this.persons.get(personID);
	}

	public PseudoPerson getPseudoPerson(String pseudopersonID) {
		return this.pseudopersons.get(pseudopersonID);
	}

	public Hashtable<String, Person> getPersons() {
		return persons;
	}

	public Hashtable<String, PseudoPerson> getPseudoPersons() {
		return pseudopersons;
	}

	public boolean hasAncestor(Person p) {
		if (p != null) {
			if (p.getDadID().equals("0") && p.getMomID().equals("0")) {
				return false;
			} else {
				if (persons.containsKey(p.getDadID()) || persons.containsKey(p.getMomID())) {
					return true;
				} else {
					return (hasAncestor(p.getDadID()) || hasAncestor(p.getMomID()));
				}
			}
		} else {
			return false;
		}	
	}

	public boolean hasAncestor(String id) {
		Person per = (Person) persons.get(id);
		if (per != null) {
			if (per.getDadID().equals("0") && per.getMomID().equals("0")) {
				return false;
			} else {
				if (persons.containsKey(per.getDadID()) || persons.containsKey(per.getMomID())) {
					return true;
				} else {
					return (hasAncestor(per.getDadID()) || hasAncestor(per.getMomID()));
				}
			}
		} else {
			return false;
		}
	}

	public void removePerson(String id) {
		persons.remove(id);
	}

	/**
	 * sets the family name (familyName)
	 * 
	 * @param familyName
	 */
	public void setFamilyName(String familyStructName) {
		this.familyStructName = familyStructName;
	}

	public void setNumMarkers(int nMarkers) {
		numMarkers = nMarkers;
	}

	/**
	 * returns the persons Hashtable (containing persons)
	 * 
	 * @return persons Hashtable
	 */
	/**
	 * sets the persons hashtable
	 * 
	 * @param persons
	 */
	public void setPersons(Hashtable<String, Person> persons) {
		this.persons = persons;
	}

	public void setMendErrs(int mendErrors) {
		this.mendErrors = mendErrors;
	}
}