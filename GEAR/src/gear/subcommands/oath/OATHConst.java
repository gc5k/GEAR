package gear.subcommands.oath;

public class OATHConst {
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

	public static final String RAF = "RAF";
	public static int raf = 5;

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

	public static boolean isNASNP(String snp) {
		for (String naStr : naStrs) {
			if (naStr.equalsIgnoreCase(snp)) {
				return true;
			}
		}
		return false;
	}

	public static boolean isNA(String s) {
		for (String naStr : naStrs) {
			if (naStr.equalsIgnoreCase(s)) {
				return true;
			}
		}
		return false;
	}

	private static final String[] naStrs = { "-9", "na", "-Inf", "Inf", ".", "-" };
}
