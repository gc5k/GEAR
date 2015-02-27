package gear.gwassummary;

public class MetaStat
{
	public MetaStat(String snp, float effect, float se, double p, char a1, boolean logit)
	{
		this.snp = snp;
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

	public void setChr(int chr)
	{
		this.chr = chr;
	}
	
	public void setBP(long bp)
	{
		this.bp = bp;
	}
	
	public void setP(double p)
	{
		this.p = p;
	}

	public void setA2(char a2)
	{
		if (a2 >=97 && a2 <= 122)
		{
			this.a2 = (char) (a2 - 32);
		}
		else
		{
			this.a2 = a2;
		}
	}
	
	public String getSNP()
	{
		return snp;
	}

	public int getChr()
	{
		return chr;
	}
	
	public long getBP()
	{
		return bp;
	}
	
	public float getEffect()
	{
		return (float) (logit? Math.log(effect) : effect);
	}

	public float getSE()
	{
		return se;
	}

	public double getP()
	{
		return p;
	}

	public char getA1()
	{
		return a1;
	}

	public char getA2()
	{
		return a2;
	}

	private String snp;
	private int chr = -1;
	private long bp = -1;
	private float se;	
	private float effect;
	private double p;
	private char a1;
	private char a2;
	private boolean logit;
}
