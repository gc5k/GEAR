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
	public long getTests()
	{
		return tests;
	}
	public void setTests(long tests)
	{
		this.tests = tests;
	}
	public double getB()
	{
		return b;
	}
	public void setRegression(double b)
	{
		this.regFlag = true;
		this.chisqFlag = false;
		this.b = b;
	}
	public boolean getRegressionFlag()
	{
		return regFlag;
	}
	public void setChisq()
	{
		this.chisqFlag = true;
		this.regFlag = false;
	}
	public boolean getChisqFlag()
	{
		return chisqFlag;
	}
	public void setMissingRate(double m)
	{
		missingRate = m;
	}
	public double getMissingRate()
	{
		return missingRate;
	}
	private double alpha = 0.05;
	private double beta = 0.05;
	private long tests;
	private double b = 0.95;
	private boolean regFlag = false;
	private boolean chisqFlag = true;
	private double missingRate = 0.01;
}
