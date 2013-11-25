package gear.subcommands.profile;

public class QRange
{
	private String name;
	private float lowerBound;
	private float upperBound;
	
	public QRange(String name, float lowerBound, float upperBound)
	{
		this.name = name;
		this.lowerBound = lowerBound;
		this.upperBound = upperBound;
	}
	
	public String getName()
	{
		return name;
	}
	
	public float getLowerBound()
	{
		return lowerBound;
	}
	
	public float getUpperBound()
	{
		return upperBound;
	}

	public boolean isInRange(float qScore)
	{
		return lowerBound <= qScore && qScore <= upperBound;
	}
}
