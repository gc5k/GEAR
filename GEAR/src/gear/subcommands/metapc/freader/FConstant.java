package gear.subcommands.metapc.freader;

public class FConstant
{
	public static final int SNP = 0;
	public static final int CHR = 1;
	public static final int A1 = 2;
	public static final int A2 = 3;
	public static final int Fvalue = 4;
	public static final int BP = 5;
	public static final int N = 6;
	public static final int num_key=7;

	public static boolean isNASNP(String snp)
	{
		String[] naSNP = initNASNP();
		for (String naStr : naSNP)
		{
			if (naStr.equalsIgnoreCase(snp))
			{
				return true;
			}
		}
		return false;
	}

	private static String[] initNASNP()
	{
		return new String[] { ".", "-", "na", "NA", "N.A.", "n.a." };
	}

}
