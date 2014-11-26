package gear.util;

import gear.ConstValues;

public class SNPMatch
{
	public static boolean isBiallelic(char a1, char a2)
	{
		return a1 != ConstValues.MISSING_ALLELE_CHAR && a2 != ConstValues.MISSING_ALLELE_CHAR;
	}

	public static boolean isAllelesMatchForTwoLoci(char a1, char a2, char b1, char b2)
	{
		if (a1 >= 97 && a1 <= 122)
		{
			a1 -= 32;
		}
		if (a2 >= 97 && a2 <= 122)
		{
			a2 -= 32;
		}
		if (b1 >= 97 && b1 <= 122)
		{
			b1 -= 32;
		}
		if (b2 >= 97 && b2 <= 122)
		{
			b2 -= 32;
		}

		StringBuffer g1 = new StringBuffer();
		StringBuffer g2 = new StringBuffer();

		if(a1 < a2)
		{
			g1.append(a1);
			g1.append(a2);
		}
		else
		{
			g1.append(a2);
			g1.append(a1);
		}

		if(b1 < b2)
		{
			g2.append(b1);
			g2.append(b2);
		}
		else
		{
			g2.append(b2);
			g2.append(b1);
		}

		boolean isMatch = true;
		if(g1.toString().compareTo(g2.toString()) != 0)
		{
			isMatch = false;
		}
		return isMatch;
	}

	public static boolean isAllelesFlipMatchForTwoLoci(char a1, char a2, char b1, char b2)
	{
		char fa1 = Flip(a1);
		char fa2 = Flip(a2);
		boolean isMatch = isAllelesMatchForTwoLoci(fa1, fa2, b1, b2);
		return isMatch;
	}

	
	public static String Flip(String a)
	{
		String f = "A";
		if (a.compareTo("A") == 0)
		{
			f = "T";
		} 
		else if (a.compareTo("a") == 0)
		{
			f = "t";
		}
		else if (a.compareTo("C") == 0)
		{
			f = "G";
		}
		else if (a.compareTo("c") == 0)
		{
			f = "g";
		}
		else if (a.compareTo("G") == 0)
		{
			f = "C";
		}
		else if (a.compareTo("g") == 0)
		{
			f = "c";
		}
		else if (a.compareTo("T") == 0)
		{
			f = "A";
		}
		else if (a.compareTo("t") == 0)
		{
			f = "a";
		}
		else
		{
			f = ConstValues.MISSING_ALLELE_STRING;
		}
		return f;
	}

	public static char Flip(char a)
	{
		char f = 'A';
		if (a == 'A')
		{
			f = 'T';
		}
		else if (a == 'a')
		{
			f = 't';
		}
		else if (a == 'C')
		{
			f = 'G';
		}
		else if (a == 'c')
		{
			f = 'g';
		}
		else if (a == 'G')
		{
			f = 'C';
		}
		else if (a == 'g')
		{
			f = 'c';
		}
		else if (a == 'T')
		{
			f = 'A';
		}
		else if (a == 't')
		{
			f = 'a';
		}
		else
		{
			f = ConstValues.MISSING_ALLELE_CHAR;
		}
		return f;
	}

	public static boolean isAmbiguous(String a1, String a2)
	{
		return isAmbiguous(a1.charAt(0), a2.charAt(0));
	}

	public static boolean isAmbiguous(char a1, char a2)
	{		
		if (a1 > a2)
		{
			char t = a1;
			a1 = a2;
			a2 = t;
		}

		if(a1 >= 97 && a1 <= 122)
		{
			a1 -= 32;
		}
		if(a2 >= 97 && a2 <= 122)
		{
			a2 -= 32;
		}

		if (a1 == 'A' && a2 == 'T')
		{
			return true;
		}

		if (a1 == 'C' && a2 == 'G')
		{
			return true;
		}
		return false;
	}

	public static boolean IsBiallelic(char a1, char a2, char b1, char b2)
	{
		char La1 = a1;
		char La2 = a2;
		char Lb1 = b1;
		char Lb2 = b2;

		boolean f = true;
		if ((La1 - La2) > 0)
		{
			char t = La1;
			La1 = La2;
			La2 = t;
		}
		StringBuilder a = new StringBuilder(La1);
		a.append(La2);
		String g1 = a.toString();

		if ((Lb1 - Lb2) > 0)
		{
			char t = Lb1;
			Lb1 = Lb2;
			Lb2 = t;
		}
		StringBuilder b = new StringBuilder(Lb1);
		b.append(Lb2);
		String g2 = b.toString();

		if (La1 != La2 && Lb1 != Lb2)
		{ // both bialelic
			if (g1.compareTo("AC") == 0)
			{
				if (g2.compareTo("AC") == 0 || g2.compareTo("GT") == 0)
				{
					return true;
				}
			}
			else if (g1.compareTo("AG") == 0)
			{
				if (g2.compareTo("AG") == 0 || g2.compareTo("CT") == 0)
				{
					return true;
				}
			}
			else if (g1.compareTo("AT") == 0)
			{
				if (g2.compareTo("AT") == 0)
				{
					return true;
				}
			}
			else if (g1.compareTo("CG") == 0)
			{
				if (g2.compareTo("CG") == 0)
				{
					return true;
				}
			}
			else if (g1.compareTo("CT") == 0)
			{
				if (g2.compareTo("CT") == 0 || g2.compareTo("AG") == 0)
				{
					return true;
				}
			}
			else
			{
				if (g2.compareTo("GT") == 0 || g2.compareTo("AC") == 0)
				{
					return true;
				}
			}
		} 
		else
		{
			return false;
		}
		return f;
	}
}
