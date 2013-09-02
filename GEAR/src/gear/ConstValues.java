package gear;

public class ConstValues
{
	
	public static byte PLINK_BED_BYTE1 = 0x6c;
	public static byte PLINK_BED_BYTE2 = 0x1b;
	public static byte PLINK_BED_BYTE3 = 0x1;

	public static final int FLOAT_SIZE = Float.SIZE / Byte.SIZE;
	
	public static final String WHITESPACE_DELIMITER = "\\s+";
	
	public static final int BINARY_HOMOZYGOTE_FIRST = 0x0;
	public static final int BINARY_HETEROZYGOTE = 0x1;
	public static final int BINARY_HOMOZYGOTYE_SECOND = 0x2;
	public static final int BINARY_MISSING_GENOTYPE = 0x3;
	
	public static final char MISSING_ALLELE_CHAR = '0';
	public static final String MISSING_ALLELE_STRING = "0";
	
	public static boolean isNA(String s)
	{
		for (String naStr : naStrs)
		{
			if (naStr.equals(s))
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
			return new String[] { "-9", "NA", "na", "-Inf", "Inf" };
		}
		return CmdArgs.INSTANCE.getNA().trim().split(",");
	}
	
	public static final double EPSILON = 1e-6;
	
	public static final long KILOBYTE = 1024;
	public static final long MEGABYTE = 1024 * 1024;
	public static final long GIGABYTE = 1024 * 1024 * 1024;
}
