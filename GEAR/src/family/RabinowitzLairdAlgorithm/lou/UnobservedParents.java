package family.RabinowitzLairdAlgorithm.lou;

import java.util.Iterator;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.ArrayList;

import family.RabinowitzLairdAlgorithm.AbstractGenoDistribution;
import gear.CmdArgs;
import gear.util.NewIt;

/**
 * Class extends GenoDistribution when neither parental genotype are available
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class UnobservedParents extends AbstractGenoDistribution
{

	/**
	 * Construct UnobservedParents
	 * 
	 * @param child
	 *            genotypes of observed kids
	 */
	public UnobservedParents(TreeMap<String, Integer> child)
	{
		super(child);
		countAllele(childrenGenoMap);
		genotypeParents();
	}

	public String[] getNontransmitted()
	{
		return null;
	}

	/**
	 * Get nontransmitted genotypes If transmitted genotye is missing and can
	 * not be randomly assigned, the nontransmitted genotype will be considered
	 * missing too.
	 */
	public String[] getNontransmitted(final String transmitted)
	{
		String nontran = new String(CmdArgs.INSTANCE.missingGenotype);
		String tran = new String(transmitted);

		if (transmitted.compareTo(CmdArgs.INSTANCE.missingGenotype) == 0)
		{
			tran = RandomAssign();
			if (tran.compareTo(CmdArgs.INSTANCE.missingGenotype) == 0)
			{
				String nontran_tran[] = { CmdArgs.INSTANCE.missingGenotype,
						CmdArgs.INSTANCE.missingGenotype };
				return nontran_tran;
			}
		}
		if (isParentGenotyped())
		{
			String nontran_tran[] = produceNontransmitted(tran);
			// System.out.println(nontran_tran[0]+"   "+nontran_tran[1]);
			return nontran_tran;
		} else if (childrenGenoMap.size() == 1)
		{// situation 1,2

			nontran = new String((String) childrenGenoMap.firstKey());
		} else if (childrenGenoMap.size() == 2)
		{
			if (numHomozygous(childrenGenoMap) > 0)
			{// situation 3

				nontran = (tran.compareTo(childrenGenoMap.firstKey()) != 0) ? childrenGenoMap
						.firstKey() : childrenGenoMap.lastKey();

				if (!isHeterozygous(transmitted))
				{// 3-1
					nontran = (transmitted
							.compareTo(childrenGenoMap.firstKey()) != 0) ? childrenGenoMap
							.firstKey() : childrenGenoMap.lastKey();
				} else
				{// 3-2
					double vk1 = childrenGenoMap
							.get(childrenGenoMap.firstKey()).doubleValue();
					double vk2 = childrenGenoMap.get(childrenGenoMap.lastKey())
							.doubleValue();
					double rd = rnd.nextFloat() * (vk1 + vk2);
					nontran = rd < vk1 ? childrenGenoMap.firstKey()
							: childrenGenoMap.lastKey();
				}
			} else
			{// situation 6, 7

				nontran = tran.compareTo(childrenGenoMap.firstKey()) != 0 ? childrenGenoMap
						.firstKey() : childrenGenoMap.lastKey();
			}
		} else
		{// situation 12

			ArrayList<String> Geno = NewIt.newArrayList();
			for (String g : childrenGenoMap.keySet())
			{
				Geno.add(g);
			}
			int index = (int) rnd.nextFloat() * childrenGenoMap.size();
			nontran = (String) Geno.get(index);
		}
		add(nontran);

		String nontran_tran[] = { nontran, tran };
		// System.out.println(nontran_tran[0]+"   "+nontran_tran[1]);
		return nontran_tran;
	}

	protected void genotypeParents()
	{
		ArrayList<TreeSet<String>> PG = NewIt.newArrayList();
		PG.add(new TreeSet<String>());
		PG.add(new TreeSet<String>());

		Iterator<String> it = childrenGenoMap.keySet().iterator();
		String geno;
		if (numAllele() == 4 && childrenGenoMap.size() >= 3)
		{
			// three genotype is sufficient to recreat parents'genotypes
			geno = it.next();// the first genotype

			PG.get(0).add(new String(geno.substring(0, 1)));
			PG.get(1).add(new String(geno.substring(1, 2)));

			geno = it.next();
			if (!PG.get(0).contains(geno.substring(0, 1))
					&& !PG.get(0).contains(geno.substring(1, 2))
					&& !PG.get(1).contains(geno.substring(0, 1))
					&& !PG.get(1).contains(geno.substring(1, 2)))
			{// for
				// situation
				// 14
				// go to next genotype if this genotype does not match any
				// allele being assigned in current stage.

				geno = it.next();
			}
			int index = 0;
			if (PG.get(0).contains(geno.substring(0, 1))
					|| PG.get(0).contains(geno.substring(1, 2)))
			{
				index = 1;
			}
			PG.get(index)
					.add(new String((PG.get(1 - index).contains(geno.substring(
							0, 1))) ? geno.substring(1, 2) : geno.substring(0,
							1)));
			geno = it.next();
			PG.get(1 - index)
					.add(new String((PG.get(index).contains(geno
							.substring(0, 1))) ? geno.substring(1, 2) : geno
							.substring(0, 1)));
		}
		if (numAllele() == 3)
		{
			if (numHomozygous(childrenGenoMap) > 0)
			{
				for (; it.hasNext();)
				{
					geno = (String) it.next();
					if (!isHeterozygous(geno))
					{
						PG.get(0).add(geno.substring(0, 1));
						PG.get(1).add(geno.substring(0, 1));
					}
				}
				// System.out.println((String) PG.get(0).first());
				String HomoAllele = (String) PG.get(0).first();
				int index = 0;
				for (String al : getAlleleSet())
				{
					if (al.contains(HomoAllele.substring(0, 1)))
					{
						continue;
					}
					PG.get(index++).add(al);
				}
			}
		}
		if (numAllele() == 2)
		{
			if (numHomozygous(childrenGenoMap) == 2)
			{
				for (String al : getAlleleSet())
				{
					PG.get(0).add(al);
					PG.get(1).add(al);
				}
			}
		}
		if (PG.get(0).size() == 2 && PG.get(1).size() == 2)
		{
			parentGeno.add(PG.get(0).first() + PG.get(0).last());
			parentGeno.add(PG.get(1).first() + PG.get(1).last());
			// System.out.println( "P1"+(String) parentGeno.firstElement() );
			// System.out.println( "P2"+(String) parentGeno.lastElement() );
		}
	}

	public String[] produceNontransmitted(String transmitted)
	{
		char nontran[] = new char[2];
		String p1 = parentGeno.get(0);
		String p2 = parentGeno.get(1);
		char PG[][] = { { p1.charAt(0), p1.charAt(1) },
				{ p2.charAt(0), p2.charAt(1) } };
		char allele[] = new char[2];

		allele[0] = transmitted.charAt(0);
		allele[1] = transmitted.charAt(1);

		int L0[][] = new int[2][2];
		int L1[][] = new int[2][2];
		int K0 = 0;
		int K1 = 0;
		int C0 = 0;
		int C1 = 0;
		for (int i = 0; i < PG.length; i++)
		{
			for (int j = 0; j < PG[i].length; j++)
			{
				if (allele[0] == PG[i][j])
				{
					L0[i][j] = 1;
					C0++;
				}
				if (allele[1] == PG[i][j])
				{
					L1[i][j] = 1;
					C1++;
				}
			}
		}

		if (((L0[0][0] + L0[0][1]) > 0) && ((L0[1][0] + L0[1][1]) > 0))
		{
			// found the allele in p1 and p2

			K0 = 2;
		} else if ((L0[0][0] + L0[0][1]) > 0)
		{// found the allele in p1;

			K0 = 0;
		} else
		{// found the allele in p2;

			K0 = 1;
		}

		if (((L1[0][0] + L1[0][1]) > 0) && ((L1[1][0] + L1[1][1]) > 0))
		{
			// found the allele in p1 and p2

			K1 = 2;
		} else if ((L1[0][0] + L1[0][1]) > 0)
		{// found the allele in p1;

			K1 = 0;
		} else
		{// found the allele in p2;

			K1 = 1;
		}

		/**
		 * All possible combination here are
		 * {(C0,C1)|(1,1),(1,2),(1,3),(2,2),(3,3)}
		 */
		if (C0 == 1 || C1 == 1)
		{// one allele meets one same allele in parent
			// genotype, the other more than one;

			if (C0 == 1)
			{
				K1 = 1 - K0;
				nontran[0] = (L0[K0][0] == 0) ? PG[K0][0] : PG[K0][1];
				nontran[1] = (L1[K1][0] == 0) ? PG[K1][0] : PG[K1][1];
			} else
			{
				K0 = 1 - K1;
				nontran[0] = (L0[K0][0] == 0) ? PG[K0][0] : PG[K0][1];
				nontran[1] = (L1[K1][0] == 0) ? PG[K1][0] : PG[K1][1];
			}
		} else
		{// c > 1 && c2 >1

			nontran[0] = (L0[0][0] == 0) ? PG[0][0] : PG[0][1];
			nontran[1] = (L1[1][0] == 0) ? PG[1][0] : PG[1][1];
		}

		if (nontran[1] < nontran[0])
		{
			char temp = nontran[0];
			nontran[0] = nontran[1];
			nontran[1] = temp;
		}

		String nontran_tran[] = { new String(nontran), transmitted };
		return nontran_tran;
	}
}