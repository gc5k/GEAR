package gear.family.pedigree.file;

import java.text.DecimalFormat;

import gear.ConstValues;
import gear.util.Logger;

/**
 * 
 * @author Guo-Bo Chen chenguobo@gmail.com Thanks Jelai Wang.
 */
public class SNP implements Comparable<SNP>
{
	private String chr = "";
	private String name = "";
	private float distance = 0;
	private int position = 0;
	private String minor;
	private char[] alleles;
	private double[] freq;

	public SNP(String n)
	{
		name = n;
	}

	public SNP(String c, String n, float d, int p)
	{
		chr = c;
		name = n;
		distance = d;
		position = p;
	}
	
	public SNP(String name, char allele1, char allele2)
	{
		this.name = name;
		this.alleles = new char[] { allele1, allele2 };
	}

	public SNP(String c, String n, float d, int p, char allele1, char allele2)
	{
		chr = c;
		name = n;
		distance = d;
		position = p;
		alleles = new char[] { allele1, allele2 };
	}

	/**
	 * Returns the name of the SNP.
	 */
	public String getName()
	{
		return name;
	}

	/**
	 * Returns the chromosome where the SNP is located.
	 */
	public String getChromosome()
	{
		return chr;
	}

	/**
	 * Returns the physical position on the chromosome where the SNP is located.
	 */
	public int getPosition()
	{
		return position;
	}

	/**
	 * Returns the genetic distance on the chromosome where the SNP is located.
	 */
	public float getDistance()
	{
		return distance;
	}

	public void setAllele(char[] a, short[] freq)
	{
		if (a.length >= 3)
		{
			Logger.printUserError("There're more than 2 alleles of " + name
					+ ".");
			System.exit(1);
		}
		alleles = new char[2];
		System.arraycopy(a, 0, alleles, 0, 2);

		this.freq = new double[freq.length];
		System.arraycopy(freq, 0, this.freq, 0, freq.length);
	}

	public void setAllele(double[] freq)
	{
		this.freq = new double[freq.length];
		System.arraycopy(freq, 0, this.freq, 0, freq.length);
	}

	public void setAllelePolymorphism(char[] a)
	{
		alleles = new char[2];
		System.arraycopy(a, 0, alleles, 0, 2);
	}

	public String getPolymorphism(String g)
	{
		StringBuffer sb = new StringBuffer();
		switch (Integer.parseInt(g))
		{
		case 0:
			sb.append(alleles[0]);
			sb.append(alleles[0]);
			break;
		case 1:
			sb.append(alleles[0]);
			sb.append(alleles[1]);
			break;
		case 2:
			sb.append(alleles[1]);
			sb.append(alleles[0]);
		case 3:
			sb.append(alleles[1]);
			sb.append(alleles[1]);
		}
		return sb.toString();
	}

	public char[] getSNP()
	{
		return alleles;
	}

	public String toString()
	{
		DecimalFormat fmt = new DecimalFormat("#.###E0");

		StringBuffer sb = new StringBuffer();
		sb.append(name + ", ");
		sb.append(chr + ", ");
		sb.append(position + ", ");
		double d0 = 0;
		double d1 = 0;
		if (freq.length > 1)
		{
			d0 = freq[0];
			d1 = freq[1];
		} else
		{
			d0 = freq[0];
		}
		if (d0 < d1)
		{
			sb.append(alleles[0] + ", " + fmt.format(d0 / (d0 + d1)) + ", ");
		} else
		{
			sb.append(alleles[1] + ", " + fmt.format(d1 / (d0 + d1)) + ", ");
		}
		if (minor != null)
		{
			sb.append(minor);
		}
		return sb.toString();
	}

	public boolean isMonopolic()
	{
		boolean flag = false;
		if (alleles[0] == ConstValues.MISSING_ALLELE_CHAR || alleles[1] == ConstValues.MISSING_ALLELE_CHAR)
		{
			flag = true;
		}
		return flag;
	}

	public char getFirstAllele()
	{
		return alleles[0];
	}

	public char getSecAllele()
	{
		return alleles[1];
	}

	@Override
	public int compareTo(SNP snp)
	{
		if (snp.chr.compareTo(this.chr) == 0)
		{
			return this.position - snp.getPosition() ;
		} 
		else if (this.chr.compareTo(snp.chr) < 0)
		{
			return Integer.MIN_VALUE;
		}
		else
		{
			return Integer.MAX_VALUE;
		}
	}
}
