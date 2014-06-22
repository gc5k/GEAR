package gear.subcommands.lambdaD;

public class MetaStat
{
	protected MetaStat(String snp, int chr, long bp, float effect, float se, double p, char a1, boolean logit)
	{
		this.snp = snp;
		this.chr = chr;
		this.bp = bp;
		this.effect = effect;
		this.se = se;
		this.p = p;
		if (a1 >=97 && a1 <= 122)
		{
			this.a1 = (char) (a1 - 32);
		}
		else
		{
			this.a1 = a1;
		}
		this.logit = logit;
	}

	protected MetaStat(String snp, int chr, long bp, float effect, float se, double p, char a1, char a2, boolean logit)
	{
		this.snp = snp;
		this.chr = chr;
		this.bp = bp;
		this.effect = effect;
		this.se = se;
		this.p = p;
		if (a1 >=97 && a1 <= 122)
		{
			this.a1 = (char) (a1 - 32);			
		}
		else
		{
			this.a1 = a1;
		}
		
		if (a2 >=97 && a1 <= 122)
		{
			this.a2 = (char) (a2 - 32);
		}
		else
		{
			this.a2 = a2;
		}
		this.logit = logit;
	}

	protected String getSNP()
	{
		return snp;
	}

	protected int getChr()
	{
		return chr;
	}
	
	protected long getBP()
	{
		return bp;
	}
	
	protected float getEffect()
	{
		return (float) (logit? Math.log(effect) : effect);
	}

	protected float getSE()
	{
		return se;
	}

	protected double getP()
	{
		return p;
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
	private int chr;
	private long bp;
	private float se;	
	private float effect;
	private double p;
	private char a1;
	private char a2;
	private boolean logit;
}
