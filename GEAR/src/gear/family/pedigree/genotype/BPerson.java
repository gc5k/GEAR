package gear.family.pedigree.genotype;

import gear.CmdArgs;
import gear.family.plink.PLINKBinaryParser;

/**
 * stores the genotypes of each individual. this class is not thread safe
 * (untested)
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public class BPerson
{

	protected String familyID;
	protected String personID;
	protected String momID;
	protected String dadID;
	protected int gender;
	protected String affectedStatus;
	protected int numMarkers;
	protected int genoLen;
	protected int[] alleles;
	protected final int intL = 16;
	public final static int shift = 4;
	protected final int mask = 15;
	public static int MissingGenotypeCode = 3;
	public static int MissingAlleleCode = 2;

	public BPerson(int numMarkers)
	{
		this.numMarkers = numMarkers;
		if (numMarkers % intL == 0)
		{
			genoLen = numMarkers / intL;
		} else
		{
			genoLen = numMarkers / intL + 1;
		}
		alleles = new int[genoLen];
	}

	public BPerson(BPerson p)
	{
		familyID = p.getFamilyID();
		personID = p.getPersonID() + "ajhg2008";
		momID = p.getMomID();
		dadID = p.getDadID();
		gender = p.getGender();
		affectedStatus = p.getAffectedStatus();
		this.numMarkers = p.getNumMarkers();
		if (numMarkers % intL == 0)
		{
			genoLen = numMarkers / intL;
		} else
		{
			genoLen = numMarkers / intL + 1;
		}
		alleles = new int[genoLen];
	}

	/**
	 * gets the family ID
	 * 
	 * @return The familyID for this person
	 */
	public String getFamilyID()
	{
		return familyID;
	}

	/**
	 * sets the family ID
	 * 
	 * @param familyID
	 */
	public void setFamilyID(String familyID)
	{
		this.familyID = familyID;
	}

	/**
	 * gets the Person ID
	 * 
	 * @return The personID for this person
	 */
	public String getPersonID()
	{
		return personID;
	}

	/**
	 * sets the person ID
	 * 
	 * @param personID
	 */
	public void setPersonID(String personID)
	{
		this.personID = personID;
	}

	/**
	 * gets the momID for this person
	 * 
	 * @return momID
	 */
	public String getMomID()
	{
		return momID;
	}

	/**
	 * sets the momid
	 * 
	 * @param momID
	 */
	public void setMomID(String momID)
	{
		this.momID = momID;
	}

	/**
	 * gets the dad ID for this person
	 * 
	 * @return dadID
	 */
	public String getDadID()
	{
		return dadID;
	}

	/**
	 * sets the dadID
	 * 
	 * @param dadID
	 */
	public void setDadID(String dadID)
	{
		this.dadID = dadID;
	}

	/**
	 * gets the gender for this person
	 * 
	 * @return gender
	 */
	public int getGender()
	{
		return gender;
	}

	/**
	 * sets the gender
	 * 
	 * @param gender
	 */
	public void setGender(int gender)
	{
		this.gender = gender;
	}

	/**
	 * gets the affected status for this person
	 * 
	 * @return affectedStatus
	 */
	public String getAffectedStatus()
	{
		return affectedStatus;
	}

	/**
	 * sets the affected status
	 * 
	 * @param affectedStatus
	 */
	public void setAffectedStatus(String affectedStatus)
	{
		this.affectedStatus = affectedStatus;
	}

	/**
	 * returns the number of markers for this person
	 * 
	 * @return integer count of markers
	 */
	public int getNumMarkers()
	{
		return numMarkers;
	}

	public void addMarker(boolean flag, int a1, int a2, int i)
	{
		int posByte = i >> shift;
		int posBit = (i & 0xf) << 1;
		if (flag)
		{
			if (a2 == a1)
			{// add 00 or 11
				int c = (a1 + a2) << posBit;
				alleles[posByte] |= c;
			} else
			{// add 10
				int c = 1 << posBit;
				alleles[posByte] |= c;
			}
		} else
		{// add 01
			alleles[posByte] |= (3 << posBit);
		}
	}

	public void addByteGenotype(int g, int posByte, int posBit)
	{
		alleles[posByte] |= g << posBit;
	}

	public void addAllMarker(byte[] genoBytes)
	{
		for (int genoByteIdx = 0; genoByteIdx < genoBytes.length; ++genoByteIdx)
		{
			int alleleIntIdx = genoByteIdx >> 2; // one int consists of 4 bytes
			int bitPosInIntOfThisByte = (genoByteIdx & 0x3) << 3;

			// One byte can store 4 genotypes
			for (int genoIdxInByte = 0; genoIdxInByte < 4; ++genoIdxInByte)
			{
				int plinkGeno = (genoBytes[genoByteIdx] >> (genoIdxInByte << 1)) & 0x3;
				int gearGeno = PLINKBinaryParser
						.convertToGearGenotype(plinkGeno);
				alleles[alleleIntIdx] |= gearGeno << bitPosInIntOfThisByte << (genoIdxInByte << 1);
			}
		}
	}

	public String getGenotypeScoreString(int i)
	{
		int posByte = i >> shift;
		int posBit = (i & 0xf) << 1;
		int g = (alleles[posByte] >> (posBit)) & 3;
		if (g == 1)
		{// 01
			return CmdArgs.INSTANCE.missingGenotype;
		} else
		{
			if (g == 2)
			{
				return Integer.toString(1);
			} else
			{
				return Integer.toString(g);
			}
		}
	}

	public int getAlleleArrayLength()
	{
		return genoLen;
	}

	public int[] getAlleleArray()
	{
		return alleles;
	}

	public String getBiAlleleGenotypeString(int i)
	{
		int posByte = i >> shift;
		int posBit = (i & 0xf) << 1;
		int g = (alleles[posByte] >> posBit) & 3;
		if (g == 3)
		{
			return CmdArgs.INSTANCE.missingGenotype;
		} else
		{
			StringBuffer sb = new StringBuffer();
			sb.append((alleles[posByte] >> (posBit + 1)) & 1);
			sb.append(alleles[posByte] >> posBit & 1);
			return sb.toString();
		}
	}

	public int getGenotypeScore(int i)
	{
		int posByte = i >> shift;
		int posBit = (i & 0xf) << 1;
		int g = (alleles[posByte] >> (posBit)) & 3;
		return g;
		// if (g == 1) {// 01
		// return 2;
		// } else {
		// if (g == 2) {
		// return 1;
		// } else {
		// return g;
		// }
		// }
	}

	public byte getOriginalGenotypeScore(int i)
	{// this only for write back to bed file purpose
		int posByte = i >> shift;
		int posBit = (i & 0xf) << 1;
		int g = (alleles[posByte] >> (posBit)) & 3;
		switch (g)
		{
		case 0:
			g = 0;
			break;
		case 1:
			g = 2;
			break;
		case 2:
			g = 3;
			break;
		default:
			g = 1;
			break; // missing
		}
		return (byte) g;
	}

	public byte getOriginalGenotypeScore(int posByte, int posBit)
	{// this only for write back to bed file purpose
		int g = (alleles[posByte] >> (posBit)) & 3;
		switch (g)
		{
		case 0:
			g = 0;
			break;
		case 1:
			g = 2;
			break;
		case 2:
			g = 3;
			break;
		default:
			g = 1;
			break; // missing
		}
		return (byte) g;
	}

	public void setNonTransmittedGenotype(int index, String geno)
	{
		int a = Integer.parseInt(geno.substring(0, 1));
		int b = Integer.parseInt(geno.substring(1, 2));
		boolean flag = geno.compareTo(CmdArgs.INSTANCE.missingGenotype) == 0 ? false
				: true;
		addMarker(flag, a, b, index);
	}
}
