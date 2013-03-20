package util.structure;

public class Predictor2 {
	private String delim= "\\s+";
	private String SNP;
	private char A1;
	private char A2;
	private double maf;
	private double NInd;
	private String[] field;

	public Predictor2(String line, int length, int lineIdx) {
		String[] s = line.split(delim);
		if(s.length == length) {
			SNP = s[0];
			A1 = s[1].charAt(0);
			A2 = s[2].charAt(0);
			maf = Double.parseDouble(s[3]);
			NInd = Double.parseDouble(s[4]);
			field = new String[length-5];
			System.arraycopy(s, 5, field, 0, field.length);
		} else {
			System.err.println("predictor informtion incorrect at line : " + (lineIdx));
			System.exit(0);
		}
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
	
	public double getNInd() {
		return NInd;
	}

	public String getField(int idx) {
		return field[idx];
	}
	
}
