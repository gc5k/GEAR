package gear.subcommands.lambdaD;

import java.util.ArrayList;

public class GMRes implements Comparable<GMRes>
{
	private String snp;
	private int chr;
	private long bp;
	private char A1;
	private char A2;
	private int N;
	private double z;
	private double p;
	private String direct;
	private double b;
	private double se;
	private ArrayList<Integer> Int;

	public void GMRes(String snp, int chr, long bp, char A1, char A2, int N, double b, double se, ArrayList<Integer> Int)
	{
		this.snp = snp;
		this.chr = chr;
		this.bp = bp;
		this.A1 = A1;
		this.A2 = A2;
		this.N = N;
		this.b = b;
		this.se = se;
		this.Int = Int;
		z = b/se;
	}

	public int getChr()
	{
		return chr;
	}

	public long getBp()
	{
		return bp;
	}

	public String toString()
	{
		StringBuffer sb = new StringBuffer();
		sb.append(snp + " " + A1 + " " + A2 + " " + N + " " + b + " " + se + " " + Int.get(Int.size()-1));
		return sb.toString();
	}

	@Override
	public int compareTo(GMRes o)
	{
		int c = this.chr - o.getChr();
		if (c != 0)
		{
			return c;
		}
		return (new Long (this.bp - o.getBp())).intValue();
	}
}
