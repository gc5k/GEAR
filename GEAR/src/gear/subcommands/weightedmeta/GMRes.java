package gear.subcommands.weightedmeta;

public class GMRes implements Comparable<GMRes>
{
	private String snp;
	private int chr;
	private long bp;
	private char A1;
	private char A2;
	private int cohort;
	private int N;
	private double z;
	private double p;
	private String direct;
	private double b;
	private double se;

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
	
	public void SetA2(char A2)
	{
		this.A2 = A2;
	}
	
	public void SetN(int N)
	{
		this.N = N;
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

	public String printTitle()
	{
		StringBuffer sb = new StringBuffer();
		sb.append("snp chr bp A1 A2 cohort b se z p direct");
		return sb.toString();
	}

	public String toString()
	{
		StringBuffer sb = new StringBuffer();
		sb.append(snp + " " + chr + " " + bp + " " + A1 + " " + A2 + " " + cohort + " " + b + " " + se + " " + " " + z + " " + p + " " + direct);
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
