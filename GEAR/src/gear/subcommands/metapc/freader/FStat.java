package gear.subcommands.metapc.freader;

public class FStat
{
	public FStat(String snp, float effect, char A1, char A2)
	{
		this.snp = snp;
		this.effect = effect;
		if (A1 >=97 && A1 <= 122)
		{
			this.A1 = (char) (A1 - 32);
		}
		else
		{
			this.A1 = A1;
		}

		if (A2 >=97 && A2 <= 122)
		{
			this.A2 = (char) (A2 - 32);
		}
		else
		{
			this.A2 = A2;
		}
	}

	public void setChr(int chr)
	{
		this.chr = chr;
	}

	public void setBP(long bp)
	{
		this.bp = bp;
	}

	public void setN(float N)
	{
		this.n = N;
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
		return effect;
	}

	public char getA1()
	{
		return A1;
	}

	public char getA2()
	{
		return A2;
	}
	
	public float getN()
	{
		return n;
	}

	private String snp;
	private int chr = -1;
	private long bp = -1;
	private float effect;
	private char A1;
	private char A2;
	private float n;
}
