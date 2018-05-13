package gear.pscontrol.hierarchy;

import java.util.ArrayList;
import java.util.Random;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.lang3.ArrayUtils;

import gear.CmdArgs;
import gear.data.Person;
import gear.data.Family;
import gear.data.UniqueRecordSet;
import gear.family.RabinowitzLairdAlgorithm.AbstractGenoDistribution;
import gear.family.RabinowitzLairdAlgorithm.lou.HeterozygousParent;
import gear.family.RabinowitzLairdAlgorithm.lou.HomozygousParent;
import gear.family.RabinowitzLairdAlgorithm.lou.ObservedParents;
import gear.family.RabinowitzLairdAlgorithm.lou.UnobservedParents;
import gear.family.pedigree.PersonIndex;
import gear.family.pedigree.file.MapFile;
import gear.family.pedigree.file.PedigreeFile;
import gear.family.pedigree.genotype.GenoSet;
import gear.util.NewIt;

public class AJHG2008 extends ChenBase
{
	private long seed = 2010;

	public AJHG2008(PedigreeFile ped, MapFile map)
	{
		super(ped, map);
	}

	public void RevvingUp(ArrayList<PersonIndex> pt)
	{
		UniqueRecordSet<Family> families = PedData.getFamilies();
		ArrayList<Family> fam_has_sib = new ArrayList<Family>();
		PersonTable.ensureCapacity(qualified_Sib);

		ArrayList<PersonIndex> s_P = NewIt.newArrayList();

		ArrayList<Integer> SibIdx = NewIt.newArrayList();

		for (int familyIdx = 0; familyIdx < families.size(); ++familyIdx)
		{
			Family family = families.get(familyIdx);
			int si = 0;
			ArrayList<Person> plist = NewIt.newArrayList();
			for (int personIdx = 0; personIdx < family.size(); ++personIdx)
			{
				Person person = family.getPerson(personIdx);

				if (family.hasAncestor(person))
				{
					s_P.add(new PersonIndex(person, false, false));
					Person pseudoper = new Person(person);
					plist.add(pseudoper);
					s_P.add(new PersonIndex(pseudoper, true, false));
					si++;
					// s++;
				}

			}

			if (si != 0)
			{
				SibIdx.add(new Integer(si));
				for (int i = 0; i < plist.size(); i++)
				{
					family.addPerson(plist.get(i));
				}
				fam_has_sib.add(family);
			}
			// c++;
		}
		PersonTable.addAll(s_P);

		numSib = ArrayUtils.toPrimitive(SibIdx.toArray(new Integer[0]));

		AbstractGenoDistribution.rnd = new Random(seed);
		NonTransmitted(fam_has_sib);
	}

	public void setSeed(long sd)
	{
		seed = sd;
	}

	private void NonTransmitted(ArrayList<Family> families)
	{
		Person pseudoper;

		for (int markerIdx = 0; markerIdx < MapData.getMarkerNumber(); ++markerIdx)
		{
			for (Family family : families)
			{
				GenoSet genoset = GenotypeSummary(family, markerIdx);
				int numGenotypedParents = genoset.getNumTypedParents();
				AbstractGenoDistribution gDis;
				TreeSet<String> aSet = new TreeSet<String>();
				if (numGenotypedParents == 2)
				{
					family.countAllele(genoset.getchildrenGenoMap(), aSet);
					family.countAllele(genoset.getparentsGenoMap(), aSet);
					gDis = new ObservedParents(genoset.getchildrenGenoMap(),
							genoset.getparentsGenoMap());
				} 
				else
				{
					family.countAllele(genoset.getchildrenGenoMap(), aSet);
					family.countAllele(genoset.getparentsGenoMap(), aSet);
					if (numGenotypedParents == 1)
					{
						String PG = genoset.getparentsGenoMap().firstKey();
						if (!AbstractGenoDistribution.isHeterozygous(PG))
						{
							// table 1
							gDis = new HomozygousParent(
									genoset.getchildrenGenoMap(),
									genoset.getparentsGenoMap());
						} 
						else
						{// table 2
							gDis = new HeterozygousParent(
									genoset.getchildrenGenoMap(),
									genoset.getparentsGenoMap());
						}
					} 
					else
					{// table 3
						gDis = new UnobservedParents(genoset.getchildrenGenoMap());
					}
				}
				
				for (int personIdx = 0; personIdx < family.size(); ++personIdx)
				{
					Person person = family.getPerson(personIdx);
					if (person.getID().contains("ajhg2008") || !family.hasAncestor(person))
					{
						continue;
					}
					pseudoper = family.getPerson(person + "ajhg2008");
					String g = person.getBiAlleleGenotypeString(markerIdx);
					boolean f = g.compareTo(CmdArgs.INSTANCE.missingGenotype) != 0;
					String[] nontran_tran = new String[2];
					if (f)
					{
						nontran_tran = family.getNonTransmitted(g, gDis);
						int a1 = Integer.parseInt(nontran_tran[0].substring(0,
								1));
						int a2 = Integer.parseInt(nontran_tran[0].substring(1,
								2));
						pseudoper.addMarker(f, a1, a2, markerIdx);
					} 
					else
					{
						pseudoper.addMarker(f, 0, 0, markerIdx);
					}
				}
			}
		}
	}

	protected GenoSet GenotypeSummary(Family family, int idx)
	{

		TreeMap<String, Integer> Ps = NewIt.newTreeMap();
		TreeMap<String, Integer> Ks = NewIt.newTreeMap();
		for (int personIdx = 0; personIdx < family.size(); ++personIdx)
		{
			Person person = family.getPerson(personIdx);
			if (person.getID().contains("ajhg2008"))
				continue;
			Person per = family.getPerson(person.getID());
			String genotype = per.getBiAlleleGenotypeString(idx);
			if (family.hasAncestor(per.getPersonID()))
			{
				if (Ks.containsKey(genotype))
				{
					Integer c = ((Integer) Ks.get(genotype));
					Ks.put(genotype, ++c);
				} 
				else
				{
					Integer c = new Integer(1);
					Ks.put(new String(genotype), c);
				}
			} 
			else
			{
				if (Ps.containsKey(genotype))
				{
					Integer c = ((Integer) Ps.get(genotype));
					Ps.put(genotype, ++c);
				} 
				else
				{
					Integer c = new Integer(1);
					Ps.put(new String(genotype), c);
				}
			}
		}
		return new GenoSet(Ps, Ks, idx);
	}

}
