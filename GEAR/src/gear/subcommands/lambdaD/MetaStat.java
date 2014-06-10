package gear.subcommands.lambdaD;

public class MetaStat
{
	protected MetaStat(String snp, float effect, float se, char a1)
	{
		this.snp = snp;
		this.effect = effect;
		this.se = se;
		this.a1 = a1;
	}
	
	protected MetaStat(String snp, float effect, float se, char a1, char a2)
	{
		this.snp = snp;
		this.effect = effect;
		this.se = se;
		this.a1 = a1;
		this.a2 = a2;
	}

	protected String getSNP()
	{
		return snp;
	}
	
	protected float getEffect()
	{
		return effect;
	}
	
	protected float getSE()
	{
		return se;
	}

	protected char getA1()
	{
		return a1;
	}
	
	protected char getA2()
	{
		return a2;
	}

	private String snp;
	private float effect;
	private float se;
	private char a1;
	private char a2;
}
