package gear.util.structure;

import gear.util.Logger;

public class Predictor
{
	private String delim = "\\s+";
	private String SNP;
	private char A1;
	private String[] field;

	public Predictor(String line, int length, int lineIdx)
	{
		String[] s = line.split(delim);
		if (s.length == length)
		{
			SNP = s[0];
			A1 = s[1].charAt(0);
			field = new String[length - 2];
			System.arraycopy(s, 2, field, 0, field.length);
		} else
		{
			Logger.printUserError("Predictor informtion incorrect at line "
					+ lineIdx + ", expect " + length + "column(s), but read "
					+ s.length);
			System.exit(1);
		}
	}

	public String getSNP()
	{
		return SNP;
	}

	public char getA1()
	{
		return A1;
	}

	public String getField(int idx)
	{
		return field[idx];
	}
}
