package gear.data;

import gear.family.pedigree.genotype.BFamilyStruct;
import gear.util.Logger;

import java.util.ArrayList;
import java.util.Hashtable;

public class FamilySet
{
	public void putFamily(BFamilyStruct family)
	{
		if (familyMap.containsKey(family.getFamilyStructName()))
		{
			Logger.fatalInternalBug("Family '" + family.getFamilyStructName() + "' is added twice.");
		}
		familyMap.put(family.getFamilyStructName(), family);
		familyList.add(family);
	}
	
	public BFamilyStruct getFamily(String famID)
	{
		return familyMap.get(famID);
	}
	
	public BFamilyStruct getFamily(int index)
	{
		return familyList.get(index);
	}
	
	public int size()
	{
		return familyList.size();
	}
	
	private Hashtable<String, BFamilyStruct> familyMap = new Hashtable<String, BFamilyStruct>();
	private ArrayList<BFamilyStruct> familyList = new ArrayList<BFamilyStruct>();
}
