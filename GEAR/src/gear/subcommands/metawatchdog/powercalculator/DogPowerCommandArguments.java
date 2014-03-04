package gear.subcommands.metawatchdog.powercalculator;

import gear.subcommands.CommandArguments;

public class DogPowerCommandArguments extends CommandArguments
{
	public double getAlpha()
	{
		return alpha;
	}
	public void setAlpha(double alpha)
	{
		this.alpha = alpha;
	}
	public double getBeta()
	{
		return beta;
	}
	public void setBeta(double beta)
	{
		this.beta = beta;
	}
	public int getTests()
	{
		return tests;
	}
	public void setTests(int tests)
	{
		this.tests = tests;
	}
	public double getH2()
	{
		return h2;
	}
	public void setH2(double h2)
	{
		this.h2 = h2;
	}
	private double alpha;
	private double beta;
	private int tests;
	private double h2;
}
