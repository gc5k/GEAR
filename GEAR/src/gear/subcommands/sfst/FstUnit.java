package gear.subcommands.sfst;

import gear.gwassummary.MetaStat;

public class FstUnit implements Comparable<FstUnit>
{

	public FstUnit(MetaStat ms1, MetaStat ms2, boolean lineup, double n1, double n2)
	{
		this.ms1 = ms1;
		this.ms2 = ms2;
		this.se1 = ms1.getSE();
		this.se2 = ms2.getSE();

		b1 = ms1.getEffect();
		b2 = ms2.getEffect();

		this.n1 = n1;
		this.n2 = n2;
		if (!lineup) b2 = 1-b2;
		calFstBW();
	}

	private void calFstBW()
	{
		double f_m = n1 / (n1+n2) * b1 + n2 / (n1+n2) * b2;
		double s_m = (n1 + n2) / 2;
		fstBW = (n1/s_m * (b1-f_m)*(b1-f_m) + n2/s_m * (b2-f_m) * (b2-f_m))/(f_m * (1-f_m));
	}

	public double getFstBW()
	{
		return fstBW;
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

	public double getN1()
	{
		return n1;
	}
	
	public double getN2()
	{
		return n2;
	}

	private double fstBW;
	private double b1;
	private double b2;
	private double se1;
	private double se2;
	private double n1;
	private double n2;

	private MetaStat ms1;
	private MetaStat ms2;


	@Override
	public int compareTo(FstUnit o)
	{
		return new Double (this.fstBW).compareTo(new Double (o.getFstBW()));
	}

}
