package gear.subcommands.weightedmeta.util;

import java.text.DecimalFormat;
import java.util.ArrayList;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ChiSquaredDistributionImpl;

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
	private double Q;
	private double pQchisq;
	private double H2;
	private double I2 = 0;

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

	public void SetQ(ArrayList<Double> beta, ArrayList<Double> wt)
	{
		//see Guido Schwarzer et al [Meta-Analysis with R, Springer], page 34, 40
		
		if (wt.size() < 2) {
			H2 = 0;
			I2 = 0;
		} else {
			for (int i = 0; i < beta.size(); i++)
			{
				Q += wt.get(i) * (beta.get(i) - b) * (beta.get(i) - b);
			}
			ChiSquaredDistributionImpl chi = new ChiSquaredDistributionImpl(wt.size()-1);
			try
			{
				pQchisq = 1-chi.cumulativeProbability(Q);
			}
			catch (MathException e)
			{
				e.printStackTrace();
			}
			H2 = Q/(wt.size() - 1);
			if (Q > (wt.size() - 1))
			{
				I2 = (H2 - 1)/H2;
			}			
		}
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
		sb.append("SNP CHR BP A1 A2 COHORT B SE Z P Direction Q p[Q] H2 I2");
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
		
		if (p > 0.001)
		{
			sb.append(fmt.format(p) + " ");
		}
		else
		{
			sb.append(fmtp.format(p) + " ");
		}
		sb.append(direct + " ");

		if (Q > 0.001)
		{
			sb.append(fmt.format(Q) + " ");
		}
		else
		{
			sb.append(fmtp.format(Q) + " ");
		}

		if (pQchisq > 0.001)
		{
			sb.append(fmt.format(pQchisq) + " ");
		}
		else
		{
			sb.append(fmtp.format(pQchisq) + " ");
		}
		
		if (H2 > 0.001)
		{
			sb.append(fmt.format(H2) + " ");
		}
		else
		{
			sb.append(fmtp.format(H2) + " ");
		}

		if (I2 > 0.001)
		{
			sb.append(fmt.format(I2));
		}
		else
		{
			sb.append(fmtp.format(I2));
		}

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
