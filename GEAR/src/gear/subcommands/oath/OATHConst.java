package gear.subcommands.oath;

import gear.CmdArgs;

public class OATHConst 
{
	public static final String SNP = "SNP";
	public static int snp = 0;

	public static final String CHR = "CHR";
	public static int chr = 1;

	public static final String BP = "BP";
	public static int bp = 2;

	public static final String RefAle = "A1";
	public static int refale = 3;

	public static final String AltAle = "A2";
	public static int altale = 4;

	public static final String Freq = "FREQ";
	public static int freq = 5;

	public static final String Vg = "VG";
	public static int vg = 6;

	public static final String BETA = "BETA";
	public static int beta = 7;

	public static final String SE = "SE";
	public static int se = 8;

	public static final String CHI = "CHI";
	public static int chi = 9;

	public static final String P = "P";
	public static int p = 10;
	
	public static final int num_key = 11;
	
	public static boolean isNASNP(String snp)
	{
		String[] naSNP = initNAStrs();
		for (String naStr : naSNP)
		{
			if (naStr.equalsIgnoreCase(snp))
			{
				return true;
			}
		}
		return false;
	}
	
	public static boolean isNA(String s)
	{
		for (String naStr : naStrs)
		{
			if (naStr.equalsIgnoreCase(s))
			{
				return true;
			}
		}
		return false;
	}

	private static final String[] naStrs = initNAStrs();

	private static String[] initNAStrs()
	{
		if (CmdArgs.INSTANCE.getNA() == null)
		{
			return new String[] { "-9", "na", "-Inf", "Inf", ".", "-"};
		}
		return CmdArgs.INSTANCE.getNA().trim().split(",");
	}

}
