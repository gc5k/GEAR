package gear.data;

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

import java.util.Set;
import java.util.TreeMap;

import gear.family.RabinowitzLairdAlgorithm.AbstractGenoDistribution;

/**
 * Storing the familyName and the members of a family from a pedigree file. This
 * class is not thread safe (untested)
 */
public class Family implements UniqueRecord
{
	// save observed genotypes;
	protected UniqueRecordSet<Person> persons = new UniqueRecordSet<Person>();
	protected String id;
	
	public Family(String id)
	{
		this.id = id;
	}
	
	@Override
	public String getID()
	{
		return id;
	}

	// this two variable were added so that accomodates Lou's Test

	/**
	 * adds a member to this family (adds to persons ArrayList)
	 * 
	 * @param ind Person to add to persons ArrayList
	 */
	public void addPerson(Person per)
	{
		persons.put(per);
	}

	/**
	 * returns the number of persons of this family
	 * 
	 * @return number of persons in this family
	 */
	public int size()
	{
		return persons.size();
	}

	/**
	 * get the Person with personID
	 * 
	 * @param personID	the personID of the person we want
	 * @return the person with matching personID
	 */
	public Person getPerson(String personID)
	{
		return persons.get(personID);
	}
	
	public Person getPerson(int personIdx)
	{
		return persons.get(personIdx);
	}

	public boolean hasPerson(String personID)
	{
		return getPerson(personID) != null;
	}

	public boolean hasAncestor(Person p)
	{
		if (p != null)
		{
			if (p.getDadID().equals("0") && p.getMomID().equals("0"))
			{
				return false;
			}
			else
			{
				// FIXME: what if I have dad but no mom?
				if (persons.has(p.getDadID()) || persons.has(p.getMomID()))
				{
					return true;
				}
				else
				{
					// FIXME: this is definitely false
					return hasAncestor(p.getDadID()) || hasAncestor(p.getMomID());
				}
			}
		}
		else
		{
			return false;
		}
	}

	public boolean hasAncestor(String personID)
	{
		return hasAncestor(persons.get(personID));
	}

	public int getNumSibs()
	{
		int numSibs = 0;
		for (int ii = 0; ii < persons.size(); ++ii)
		{
			if (hasAncestor(persons.get(ii)))
			{
				++numSibs;
			}
		}
		return numSibs;
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