package gear.subcommands.lambdaD;

import gear.gwassummary.MetaStat;

public class LamUnit implements Comparable<LamUnit>
{

	public LamUnit(MetaStat ms1, MetaStat ms2, int mode, boolean lineup, double n1, double n2)
	{
		this.ms1 = ms1;
		this.ms2 = ms2;
		this.se1 = ms1.getSE();
		this.se2 = ms2.getSE();

		b1 = ms1.getEffect();
		b2 = ms2.getEffect();

		this.n1 = n1;
		this.n2 = n2;
		this.mode = mode;
		if(mode == LambdaDCommandArguments.FRQ)
		{
			if (!lineup) b2 = 1-b2;
			calFrqChi();
			calFstBW();
		}
		else if (mode == LambdaDCommandArguments.BETA)
		{
			if (!lineup) b2 *= -1;
			calBetaChi();
		}
	}

	private void calFrqChi()
	{
		calChi();
	}
	
	private void calBetaChi()
	{
		calChi();
	}

	private void calChi()
	{
		chi1 = (b1-b2) * (b1-b2) / (se1 * se1 + se2 * se2);
	}

	private void calFstBW()
	{
		double f_m = n1 / (n1+n2) * b1 + n2 / (n1+n2) * b2;
		double s_m = (n1 + n2) / 2;
		fstBW = (n1/s_m * (b1-f_m)*(b1-f_m) + n2/s_m * (b2-f_m) * (b2-f_m))/(f_m * (1-f_m));
	}

	public double getIndicateStat(int m)
	{
		if(m == LambdaDCommandArguments.FRQ)
		{
			return getFrqChi();
		}
		else if (m == LambdaDCommandArguments.BETA)
		{
			return getBetaChi();
		}
		else
		{
			return Double.NaN;
		}
	}

	public double getBetaChi()
	{
		return chi1;
	}
	
	public double getFrqChi()
	{
		return chi1;
	}

	public double getChi()
	{
		return chi1;
	}

	public double getFstBW()
	{
		return fstBW;
	}

	public double getChi1()
	{
		return chi1;
	}
	
	public MetaStat getMetaStat1()
	{
		return ms1;
	}
	
	public MetaStat getMetaStat2()
	{
		return ms2;
	}

	public double getB1()
	{
		return b1;
	}
	
	public double getB2()
	{
		return b2;
	}

	public double getSE1()
	{
		return se1;
	}
	
	public double getSE2()
	{
		return se2;
	}

	private double chi1;
	private double fstBW;
	private double b1;
	private double b2;
	private double se1;
	private double se2;
	private double n1;
	private double n2;

	private int mode = 0;

	private MetaStat ms1;
	private MetaStat ms2;

	@Override
	public int compareTo(LamUnit o)
	{
		if(mode == LambdaDCommandArguments.BETA || mode == LambdaDCommandArguments.FRQ)
		{
			return new Double (this.chi1).compareTo(new Double (o.getChi1()));
		}
		else
		{
			return new Double (this.fstBW).compareTo(new Double (o.getFstBW()));
		}
	}
}
