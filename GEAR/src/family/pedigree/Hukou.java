package family.pedigree;

public class Hukou
{

	private String FamilyID;
	private String IndividualID;
	private String FatherID;
	private String MotherID;
	private String Sex;
	private String col6;

	private boolean isAvailable;

	public Hukou(String fid, String iid, String fa, String mo, String sex,
			String c6)
	{
		FamilyID = fid;
		IndividualID = iid;
		FatherID = fa;
		MotherID = mo;
		Sex = sex;
		col6 = c6;
	}

	public Hukou(String fid, String iid)
	{
		FamilyID = fid;
		IndividualID = iid;
	}

	public String getFamilyID()
	{
		return FamilyID;
	}

	public String getIndividualID()
	{
		return IndividualID;
	}

	public String getFatherID()
	{
		return FatherID;
	}

	public String getMotherID()
	{
		return MotherID;
	}

	public String getSet()
	{
		return Sex;
	}

	public String getCol6()
	{
		return col6;
	}

	public void setAvailable(boolean a)
	{
		isAvailable = a;
	}

	public boolean isAvailable()
	{
		return isAvailable;
	}
}
