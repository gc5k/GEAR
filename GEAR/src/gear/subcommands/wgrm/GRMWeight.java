package gear.subcommands.wgrm;

public class GRMWeight
{
	protected GRMWeight(int numValues)
	{
		this.values = new Float[numValues];
	}

	protected Float getValue(int i)
	{
		return values[i];
	}

	protected void setValue(int i, Float value)
	{
		values[i] = value;
	}

	private Float[] values;
}
