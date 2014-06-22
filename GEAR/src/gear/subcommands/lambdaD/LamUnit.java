package gear.subcommands.lambdaD;

public class LamUnit
{
	public LamUnit(double chi1, MetaStat ms1, MetaStat ms2)
	{
		this.chi1 = chi1;
		this.ms1 = ms1;
		this.ms2 = ms2;
	}
	
	public double getChi1()
	{
		return chi1;
	}
	
	public MetaStat getMS1()
	{
		return ms1;
	}
	
	public MetaStat getMS2()
	{
		return ms2;
	}

	private double chi1;
	private MetaStat ms1;
	private MetaStat ms2;
}
