package gear.subcommands.at;

import gear.subcommands.CommandArguments;

public class AtCommandArguments extends CommandArguments
{

	public void setWindow(String k)
	{
		window = Double.parseDouble(k);
		if(window < 0)
		{
			windowFlag = false;
		}
		else
		{
			windowFlag = true;			
		}
	}

	public double getWindow()
	{
		return window;
	}
	
	public boolean isWindow()
	{
		return windowFlag;
	}

	public void setRsq()
	{
		rsqFlag = true;
	}
	
	public boolean isRsq()
	{
		return rsqFlag;
	}

	public void setD()
	{
		dFlag = true;
	}
	
	public boolean isD()
	{
		return dFlag;
	}
	
	public void setDPrime()
	{
		dPrimeFlag = true;
	}

	public boolean isDPrime()
	{
		return dPrimeFlag;
	}

	private boolean dFlag = false;
	private boolean dPrimeFlag = false;
	private boolean rsqFlag = false;

	private double window = -1;
	private boolean windowFlag = false;
}
