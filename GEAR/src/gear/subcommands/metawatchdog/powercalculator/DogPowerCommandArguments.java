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
	public double getB()
	{
		return b;
	}
	public void setBValue(double b)
	{
		this.b = b;
	}
	public void setRegression(boolean f)
	{
		this.regFlag = f;
	}
	public boolean getRegression()
	{
		return regFlag;
	}
	public void setChisq(boolean f)
	{
		this.chisqFlag = true;
	}
	public void setQValue(double q)
	{
		this.q = q;
	}
	public double getQValue()
	{
		return q;
	}
	public boolean getChisq()
	{
		return chisqFlag;
	}
	private double alpha;
	private double beta;
	private int tests;
	private double b;
	private boolean regFlag;
	private boolean chisqFlag;
	private double q;
}
