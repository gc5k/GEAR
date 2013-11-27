package gear.subcommands.he.assocpower;

import gear.subcommands.CommandArguments;

public class HEAssocPowerCommandArguments extends CommandArguments
{
	public int getNumberOfMarkers()
	{
		return numberOfMarkers;
	}

	public void setNumberOfMarkers(int numberOfMarkers)
	{
		this.numberOfMarkers = numberOfMarkers;
	}

	public float getHeritabilityVariance()
	{
		return heritabilityVariance;
	}

	public void setHeritabilityVariance(float heritabilityVariance)
	{
		this.heritabilityVariance = heritabilityVariance;
	}

	public float getAlpha()
	{
		return alpha;
	}

	public void setAlpha(float alpha)
	{
		this.alpha = alpha;
	}

	public float getPower()
	{
		return power;
	}

	public void setPower(float power)
	{
		this.power = power;
	}

	private int numberOfMarkers;
	private float heritabilityVariance;
	private float alpha;
	private float power;
}
