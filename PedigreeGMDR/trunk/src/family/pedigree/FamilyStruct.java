package family.pedigree;

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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import publicAccess.PublicData;
import family.RabinowitzLairdAlgorithm.AbstractGenoDistribution;
import family.RabinowitzLairdAlgorithm.HeterozygousParent;
import family.RabinowitzLairdAlgorithm.HomozygousParent;
import family.RabinowitzLairdAlgorithm.ObservedParents;
import family.RabinowitzLairdAlgorithm.Rabinowitz0;
import family.RabinowitzLairdAlgorithm.Rabinowitz1;
import family.RabinowitzLairdAlgorithm.Rabinowitz2;
import family.RabinowitzLairdAlgorithm.Rabinowitz3;
import family.RabinowitzLairdAlgorithm.UnobservedParents;

/**
 * Storing the familyName and the members of a family from a pedigree file.
 * This class is not thread safe (untested)
 * 
 * @author Julian Maller, extended by Guo-Bo Chen
 */
public class FamilyStruct {
//save observed genotypes;
    private Hashtable persons;
    private Hashtable pseudopersons;
    private String familyStructName;
    private int mendErrors;
    private int numMarkers;
    private ArrayList ObservedGenoSet;
    private ArrayList ImputedGenoSet;
    private ArrayList nontransmitted;

    /**
     * add a transmitted genoset into this family ( adds to transmitted ArrayList)
     * 
     * @param geno
     */
    public void addObservedGenoSet(GenoSet geno) {
        this.ObservedGenoSet.add(geno);
    }

    /**
     * add a genoset after imputation into this family ( adds to transmitted ArrayList)
     *
     * @param geno
     */
    public void addImputedGenoSet(GenoSet geno) {
        this.ImputedGenoSet.add(geno);
    }
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

    public void countAllele(TreeMap Geno, Set alleleSet) {
        Set GSet = Geno.keySet();
        Iterator it = GSet.iterator();
        for (; it.hasNext();) {
            String g = (String) it.next();
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
        this.persons = new Hashtable();
        this.pseudopersons = new Hashtable();
        this.ObservedGenoSet = new ArrayList();
        this.ImputedGenoSet = new ArrayList();
        this.nontransmitted = new ArrayList();
    }

    public FamilyStruct(String familyStructName) {
        this.persons = new Hashtable();
        this.pseudopersons = new Hashtable();
        this.ObservedGenoSet = new ArrayList();
        this.ImputedGenoSet = new ArrayList();
        this.nontransmitted = new ArrayList();
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
     * Table 1: Conditional distribution of non-transmitted genotypes when one homozygous parents's genotype and
     * children's genotypes are available at marker locus i Table 2: Conditional contribution of non-transmitted
     * genotypes when one heterozygous parent's genotype and children's genotype are available at marker locus i Table
     * 3: Conditional contribution of non-transmitted genotypes when only children's genotypes but no parent's genotype
     * are available at marker locus A
     * 
     * @param genotype
     *            , transmitted genotype at locus i of indiviudal j
     * @param genoset
     *            , configuration of parents' and children's genotypes.
     * @return non-transmitted genotype
     */
    public String[] getNonTransmitted(String transmitted,
            AbstractGenoDistribution genodis) {
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

    /**
     * returns a list of personIDs that are persons of this family in the form of an enumeration
     * 
     * @return enumeration memberlist
     */
    public Enumeration getPersonList() {
        return this.persons.keys();
    }

    public String[] getPersonListSorted() {
    	Enumeration perstrList = getPersonList();
    	String[] PID = new String[persons.size()];
		int ind = 0;
		while (perstrList.hasMoreElements()) {
			PID[ind++] = (String) perstrList.nextElement();
		}
		Arrays.sort(PID);
		return PID;
    }

    /**
     * get transmitted genoset
     * 
     * @return
     */
    public ArrayList getObvservedGenoSet() {
        return ObservedGenoSet;
    }

    /**
     * get transmitted genoset after imputation
     *
     * @return
     */
    public ArrayList getImputedGenoSet() {
        return ImputedGenoSet;
    }

    /**
     * get transmtted genoset at position i
     * 
     * @param i
     * @return
     */
    public GenoSet getObservedGenoSet(int i) {
        return (GenoSet) ObservedGenoSet.get(i);
    }

    /**
     * get the Person with personID
     * 
     * @param personID
     *            the personID of the person we want
     * @return the person with matching personID
     */
    public Person getPerson(String personID) {
        return (Person) this.persons.get(personID);
    }

    public PseudoPerson getPseudoPerson(String pseudopersonID) {
        return (PseudoPerson) this.pseudopersons.get(pseudopersonID);
    }

    public Hashtable getPersons() {
        return persons;
    }

    public Hashtable getPseudoPersons() {
        return pseudopersons;
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

    public boolean NontransmittedProc(ArrayList markerInfor, int[] subsetMarker)
            throws FamilyStructException {
        boolean Informative = true;
        Enumeration perList;
        Person per;
        PseudoPerson pseudoper;
        String iid;
        String nontran_tran[];
        for (int i = 0; i < subsetMarker.length; i++) {
            perList = persons.keys();
            GenoSet genoset = (GenoSet) ImputedGenoSet.get(subsetMarker[i]);
            int numGenotypedParents = genoset.getNumTypedParents();
            AbstractGenoDistribution gDis;
            TreeSet aSet = new TreeSet();
            if (genoset.getNumParents() > 2) {
                throw new FamilyStructException(
                        "Family " + familyStructName + " is not a nuclear family. It has more than 2 founders");
            }
            if (numGenotypedParents == 2) {
                countAllele(genoset.getchildrenGenoMap(), aSet);
                countAllele(genoset.getparentsGenoMap(), aSet);
                if (aSet.size() > 4) {
                    throw new FamilyStructException("Marker " + markerInfor.get(subsetMarker[i]) + " has more than 4 alleles.");
                }
                gDis = new ObservedParents(genoset.getchildrenGenoMap(),
                        genoset.getparentsGenoMap());
            } else {
                countAllele(genoset.getchildrenGenoMap(), aSet);
                countAllele(genoset.getparentsGenoMap(), aSet);
                if (numGenotypedParents == 1) {
                    String PG = (String) ((TreeMap) genoset.getparentsGenoMap()).firstKey();
                    if (!AbstractGenoDistribution.isHeterozygous(PG)) // table 1
                    {
                        if (aSet.size() > 3) {
                            throw new FamilyStructException(
                                    "Marker " + markerInfor.get(subsetMarker[i]) + " has more than 3 alleles with one parent is homozygous.");
                        }
                        gDis = new HomozygousParent(genoset.getchildrenGenoMap(), genoset.getparentsGenoMap());
                    } else // table 2
                    {
                        if (aSet.size() > 4) {
                            throw new FamilyStructException(
                                    "Marker " + markerInfor.get(subsetMarker[i]) + " has more than 4 alleles with one parent is heterozygous.");
                        }
                        gDis = new HeterozygousParent(genoset.getchildrenGenoMap(), genoset.getparentsGenoMap());
                    }
                } else // table 3
                {
                    gDis = new UnobservedParents(genoset.getchildrenGenoMap());
                }
            }
            while (perList.hasMoreElements()) {
                iid = (String) perList.nextElement();
                per = (Person) (persons.get(iid));
                pseudoper = (PseudoPerson) pseudopersons.get(iid);
                if (!hasAncestor(per.getPersonID())) {
                    continue;
                }
                nontran_tran = getNonTransmitted(per.getGenotype(subsetMarker[i]), gDis);
                pseudoper.addMarker(nontran_tran[0]);
                if (nontran_tran[0].compareTo(PublicData.MissingGenotype) == 0) {
                    Informative = false;
                    System.err.println("Individual " + per.getPersonID() + ", in family " + per.getFamilyID() + ", failed to get the nontransmitted genotype for marker " + markerInfor.get(subsetMarker[i]));
                } else {
                    if (per.getGenotype(subsetMarker[i]).compareTo(PublicData.MissingGenotype) == 0) {
                        per.setGenotype(subsetMarker[i], nontran_tran[1]);
                    }
                }
            }
        }
        return Informative;
    }

    public void removePerson(String id) {
        persons.remove(id);
    }

    public boolean RabinowitzProc(ArrayList markerInfor, int[] subsetMarker) throws FamilyStructException {
        boolean Informative = true;
        Enumeration perList;
        Person per;
        PseudoPerson pseudoper;
        String iid;
        perList = persons.keys();
        while (perList.hasMoreElements()) {
            iid = (String) perList.nextElement();
            per = (Person) (persons.get(iid));
            pseudoper = (PseudoPerson) pseudopersons.get(iid);
            if (!hasAncestor(per.getPersonID())) {
                continue;
            }
            pseudoper.pseudoGenotypeClear();
        }

        for (int i = 0; i < subsetMarker.length; i++) {
            perList = persons.keys();
            GenoSet genoset = (GenoSet) ImputedGenoSet.get(subsetMarker[i]);
            AbstractGenoDistribution gDis;
            TreeSet aSet = new TreeSet();
            if (genoset.getNumParents() > 2) {
                throw new FamilyStructException(
                        "Family " + familyStructName + " is not a nuclear family. It has more than 2 founders");
            }
            if (genoset.getNumTypedParents() == 2) {
                countAllele(genoset.getchildrenGenoMap(), aSet);
                countAllele(genoset.getparentsGenoMap(), aSet);
                if (aSet.size() > 4) {
                    throw new FamilyStructException("Marker " + markerInfor.get(subsetMarker[i]) + " has more than 4 alleles.");
                }
                gDis = new Rabinowitz0(genoset.getchildrenGenoMap(), genoset.getparentsGenoMap());
            } else {
                countAllele(genoset.getchildrenGenoMap(), aSet);
                countAllele(genoset.getparentsGenoMap(), aSet);
                if (genoset.getNumTypedParents() == 1) {
                    String PG = (String) ((TreeMap) genoset.getparentsGenoMap()).firstKey();
                    if (!AbstractGenoDistribution.isHeterozygous(PG)) {// table 1                  
                        if (aSet.size() > 3) {
                            throw new FamilyStructException(
                                    "Marker " + markerInfor.get(subsetMarker[i]) + " has more than 3 alleles with one parent is homozygous.");
                        }
                        gDis = new Rabinowitz1(genoset.getchildrenGenoMap(),
                                genoset.getparentsGenoMap());
                    } else {// table 2
                        if (aSet.size() > 4) {
                            throw new FamilyStructException(
                                    "more than 4 alleles with one parent is heterozygous.");
                        }
                        gDis = new Rabinowitz2(genoset.getchildrenGenoMap(),
                                genoset.getparentsGenoMap());
                    }
                } else {// table 3
                    gDis = new Rabinowitz3(genoset.getchildrenGenoMap());
                }
            }
            String[] controlGenotype = gDis.getNontransmitted();
            int index = 0;
            while (perList.hasMoreElements()) {
                iid = (String) perList.nextElement();
                per = (Person) (persons.get(iid));
                pseudoper = (PseudoPerson) pseudopersons.get(iid);
                if (!hasAncestor(per.getPersonID())) {
                    continue;
                }
                pseudoper.addMarker(controlGenotype[index++]);
            }
        }
        return Informative;
    }

    /**
     * sets the family name (familyName)
     * 
     * @param familyName
     */
    public void setFamilyName(String familyStructName) {
        this.familyStructName = familyStructName;
    }

    public void setNumMarkers(int numMarkers) {
        this.numMarkers = numMarkers;
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
    public void setPersons(Hashtable persons) {
        this.persons = persons;
    }

    public void setMendErrs(int mendErrors) {
        this.mendErrors = mendErrors;
    }
}