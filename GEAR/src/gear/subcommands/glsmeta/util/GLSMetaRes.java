package gear.subcommands.glsmeta.util;

import gear.util.stat.PrecisePvalue;

import org.apache.commons.math.linear.RealMatrix;

public class GLSMetaRes implements Comparable<GLSMetaRes>
{
	private String snp;
	private int chr;
	private long bp;
	private char A1;
	private char A2;
	private int cohort;
	private String direct;
	private boolean isAmbiguous;
	private RealMatrix B;
	private RealMatrix V;

	public GLSMetaRes(int cohort)
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

	public void SetBP(long bp)
	{
		this.bp = bp;
	}

	public void SetB(RealMatrix B)
	{
		this.B = B;
	}

	public void SetSE(RealMatrix V)
	{
		this.V = V;
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
		sb.append("SNP CHR BP TEST A1 A2 COHORT B SE Z P Direction");
		return sb.toString();
	}

	public String toString()
	{
		StringBuffer sb = new StringBuffer();
		for(int i = 0; i < B.getRowDimension(); i++)
		{
			double b = B.getEntry(i, 0);
			double se = Math.sqrt(V.getEntry(i, i));
			double z = b/se;
			double p = PrecisePvalue.getPvalue4Z_2Tail(z);
			if (i == 0)
			{
				sb.append(snp + " " + chr + " " + bp + " " +"SNP " + A1 + " " + A2 + " " + cohort + " " + b + " " + se + " " + " " + z + " " + p + " " + direct + "\n");				
			}
			else
			{
				sb.append(snp + " " + chr + " " + bp + " " +"COV" + i + " " + A1 + " " + A2 + " " + cohort + " " + b + " " + se + " " + " " + z + " " + p + " " + direct);				
			}
		}
		return sb.toString();
	}

	@Override
	public int compareTo(GLSMetaRes o)
	{
		int c = this.chr - o.getChr();
		if (c != 0)
		{
			return c;
		}
		return (new Long (this.bp - o.getBp())).intValue();
	}
}
