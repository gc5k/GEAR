package gear.subcommands.lambdaD;

public class MetaStat
{
	protected MetaStat(String snp, float effect, float se, char a1, boolean logit)
	{
		this.snp = snp;
		this.effect = effect;
		this.se = se;
		this.a1 = a1;
		this.logit = logit;
	}
	
	protected MetaStat(String snp, float effect, float se, char a1, char a2, boolean logit)
	{
		this.snp = snp;
		this.effect = effect;
		this.se = se;
		this.a1 = a1;
		this.a2 = a2;
		this.logit = logit;
	}

	protected String getSNP()
	{
		return snp;
	}

	protected float getEffect()
	{
		return (float) (logit? Math.log(effect) : effect);
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
	private boolean logit;
}
