package gear.subcommands.profile;

/**
 * Score values for a specific allele of a SNP 
 * 
 * @author Zhixiang
 *
 */
class Score
{
	protected Score(char allele, int numValues)
	{
		this.allele = allele;
		this.values = new Float[numValues];
	}

	protected char getAllele()
	{
		return allele;
	}

	protected Float getValue(int i)
	{
		return values[i];
	}
	
	protected void setValue(int i, Float value)
	{
		values[i] = value;
	}
	
	private char allele;
	private Float[] values;
}
