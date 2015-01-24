package gear.subcommands.randmeta;

public class BetaVec
{
	public BetaVec(double b1, double b2, double se1, double se2)
	{
		this.b1 = b1;
		this.b2 = b2;
		this.se1 = se1;
		this.se2 = se2;
		
		z1 = b1/se1;
		z2 = b2/se2;
	}

	public double getB1()
	{
		return b1;
	}
	
	public double getB2()
	{
		return b2;
	}
	
	public double getZ1()
	{
		return z1;
	}
	
	public double getZ2()
	{
		return z2;
	}

	private double b1;
	private double se1;
	private double z1;

	private double b2;
	private double se2;
	private double z2;
}
