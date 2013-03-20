package util.structure;

public class MAF {
	private String delim= "\\s+";
	private String chr;
	private String SNP;
	private char A1;
	private char A2;
	private double maf;
	private double NChr;

	public MAF(String line, int idx) {

		String[] s = line.split(delim);
		if(s.length== 6) {
			chr = s[0];
			SNP = s[1];
			A1 = s[2].charAt(0);
			A2 = s[3].charAt(0);
			maf = Double.parseDouble(s[4]);
			NChr = Double.parseDouble(s[5]);
		} else {
			System.err.println("maf informtion incorrect at line : " + (idx));
			System.exit(0);
		}
	}

	public String getChr() {
		return chr;
	}

	public String getSNP() {
		return SNP;
	}

	public char getA1() {
		return A1;
	}

	public char getA2() {
		return A2;
	}

	public double getMAF() {
		return maf;
	}

	public double getNChr() {
		return NChr;
	}
}
