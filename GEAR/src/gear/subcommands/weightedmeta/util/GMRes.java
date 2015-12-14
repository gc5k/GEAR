package gear.subcommands.weightedmeta.util;

import java.text.DecimalFormat;

public class GMRes implements Comparable<GMRes>
{
	private String snp;
	private int chr;
	private long bp;
	private char A1;
	private char A2;
	private int cohort;
	private double z;
	private double p;
	private String direct;
	private double b;
	private double se;
	private boolean isAmbiguous;

	public GMRes(int cohort)
	{
		this.cohort = cohort;
	}

	public void SetSNP(String snp)
	{
		this.snp = snp; 
	}
	
	public void SetChr(int chr)
	{
		this.chr = chr;
	}
	
	public void SetA1(char A1)
	{
		this.A1 = A1;
	}
	
	public char GetA1()
	{
		return A1;
	}
	
	public char GetA2()
	{
		return A2;
	}

	public void SetA2(char A2)
	{
		this.A2 = A2;
	}
	
	public void SetZ(double z)
	{
		this.z = z;
	}
	
	public void SetBP(long bp)
	{
		this.bp = bp;
	}

	public void SetB(double b)
	{
		this.b = b;
	}

	public void SetSE(double se)
	{
		this.se = se;
	}

	public void SetP(double p)
	{
		this.p = p;
	}

	public void SetAmbi(boolean ambi)
	{
		isAmbiguous = ambi;
	}

	public void SetDirect(String direct)
	{
		this.direct = direct;
	}

	public int getChr()
	{
		return chr;
	}

	public long getBp()
	{
		return bp;
	}

	public boolean getIsAmbiguous()
	{
		return isAmbiguous;
	}

	public String printTitle()
	{
		StringBuffer sb = new StringBuffer();
		sb.append("SNP CHR BP A1 A2 COHORT B SE Z P Direction");
		return sb.toString();
	}

	public String toString()
	{
		DecimalFormat fmt = new DecimalFormat("0.0000");
		DecimalFormat fmtp = new DecimalFormat("0.000E000");

		StringBuffer sb = new StringBuffer();
		sb.append(snp + " " + chr + " " + bp + " " + A1 + " " + A2 + " " + cohort + " ");
		if (Math.abs(b) >= 0.0001)
		{
			sb.append(fmt.format(b) + " ");
		}
		else
		{
			sb.append(fmtp.format(b) + " ");
		}
		
		if (Math.abs(se) >= 0.0001)
		{
			sb.append(fmt.format(se) + " ");
		}
		else
		{
			sb.append(fmtp.format(se) + " ");
		}
		
		if (Math.abs(z) >= 0.001)
		{
			sb.append(fmt.format(z) + " ");
		}
		else
		{
			sb.append(fmtp.format(z) + " ");
		}
		sb.append(direct);
		return sb.toString();
	}

	@Override
	public int compareTo(GMRes o)
	{
		int c = this.chr - o.getChr();
		if (c != 0)
		{
			return c;
		}
		return (new Long (this.bp - o.getBp())).intValue();
	}
}
