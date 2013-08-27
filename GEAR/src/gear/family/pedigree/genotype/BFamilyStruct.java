package gear.family.pedigree.genotype;

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

import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Set;
import java.util.TreeMap;

import gear.family.RabinowitzLairdAlgorithm.AbstractGenoDistribution;
import gear.util.NewIt;

/**
 * Storing the familyName and the members of a family from a pedigree file. This
 * class is not thread safe (untested)
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class BFamilyStruct
{
	// save observed genotypes;
	protected Hashtable<String, BPerson> persons;
	protected String familyStructName;

	// this two variable were added so that accomodates Lou's Test

	/**
	 * adds a member to this family (adds to persons ArrayList)
	 * 
	 * @param ind
	 *            Person to add to persons ArrayList
	 */
	public void addPerson(BPerson per)
	{
		this.persons.put(per.getPersonID(), per);
	}

	// public void addPseudoPerson(BPerson per) {
	// this.pseudopersons.put(per.getPersonID(), per);
	// }

	public BFamilyStruct(String familyStructName)
	{
		this.persons = NewIt.newHashtable();
		this.familyStructName = familyStructName;
	}

	/**
	 * returns the name of this family (familyName)
	 * 
	 * @return family name
	 */
	public String getFamilyStructName()
	{
		return familyStructName;
	}

	/**
	 * returns the number of persons of this family
	 * 
	 * @return number of persons in this family
	 */
	public int getNumPersons()
	{
		return this.persons.size();
	}

	/**
	 * returns a list of personIDs that are persons of this family in the form
	 * of an enumeration
	 * 
	 * @return enumeration memberlist
	 */
	public Enumeration<String> getPersonList()
	{
		return this.persons.keys();
	}

	public String[] getPersonListSorted()
	{
		Enumeration<String> perstrList = getPersonList();
		String[] PID = new String[persons.size()];
		int ind = 0;
		while (perstrList.hasMoreElements())
		{
			PID[ind++] = perstrList.nextElement();
		}
		// Arrays.sort(PID);
		return PID;
	}

	/**
	 * get the Person with personID
	 * 
	 * @param personID
	 *            the personID of the person we want
	 * @return the person with matching personID
	 */
	public BPerson getPerson(String personID)
	{
		return this.persons.get(personID);
	}

	// public BPerson getPseudoPerson(String personID) {
	// return this.pseudopersons.get(personID);
	// }

	public Hashtable<String, BPerson> getPersons()
	{
		return persons;
	}

	public boolean hasAncestor(BPerson p)
	{
		if (p != null)
		{
			if (p.getDadID().equals("0") && p.getMomID().equals("0"))
			{
				return false;
			} else
			{
				if (persons.containsKey(p.getDadID())
						|| persons.containsKey(p.getMomID()))
				{
					return true;
				} else
				{
					return (hasAncestor(p.getDadID()) || hasAncestor(p
							.getMomID()));
				}
			}
		} else
		{
			return false;
		}
	}

	public boolean hasAncestor(String id)
	{
		BPerson per = (BPerson) persons.get(id);
		if (per != null)
		{
			if (per.getDadID().equals("0") && per.getMomID().equals("0"))
			{
				return false;
			} else
			{
				if (persons.containsKey(per.getDadID())
						|| persons.containsKey(per.getMomID()))
				{
					return true;
				} else
				{
					return (hasAncestor(per.getDadID()) || hasAncestor(per
							.getMomID()));
				}
			}
		} else
		{
			return false;
		}
	}

	public int getNumSibs()
	{
		int i = 0;
		for (String pi : persons.keySet())
		{
			if (hasAncestor(pi))
			{
				i++;
			}
		}
		return i;
	}

	public void countAllele(TreeMap<String, Integer> Geno, Set<String> alleleSet)
	{

		for (String g : Geno.keySet())
		{
			Integer geno = Integer.parseInt(g);
			alleleSet.add(Integer.toString(geno & 1));
			alleleSet.add(Integer.toString(geno & 2));
		}

	}

	/**
	 * Table 1: Conditional distribution of non-transmitted genotypes when one
	 * homozygous parents's genotype and children's genotypes are available at
	 * marker locus i Table 2: Conditional contribution of non-transmitted
	 * genotypes when one heterozygous parent's genotype and children's genotype
	 * are available at marker locus i Table 3: Conditional contribution of
	 * non-transmitted genotypes when only children's genotypes but no parent's
	 * genotype are available at marker locus A
	 * 
	 * @param genotype
	 *            , transmitted genotype at locus i of indiviudal j
	 * @param genoset
	 *            , configuration of parents' and children's genotypes.
	 * @return non-transmitted genotype
	 */
	public String[] getNonTransmitted(String transmitted,
			AbstractGenoDistribution genodis)
	{
		String nontran_tran[] = genodis.getNontransmitted(transmitted);
		// System.out.println(nontran_tran[0] + "\t" + nontran_tran[1]);
		return nontran_tran;

	}
}