package gear.util.stat;

public class PrecisePvalue
{
	public static double ZcumulativeProbability(double x)
	{
		//http://en.wikipedia.org/wiki/Q-function#cite_note-6
		//IEEE Communications Letters, 2007, 11, 644
		double p1 = (1-Math.exp(-1.4*x))*Math.exp(-x*x/2)/(1.135 * Math.sqrt(2*Math.PI) * x);
		return p1*2;
	}

	public static double TwoTailZcumulativeProbability(double x)
	{
		//http://en.wikipedia.org/wiki/Q-function#cite_note-6
		//IEEE Communications Letters, 2007, 11, 644
		double p1 = ZcumulativeProbability(x);
		return p1*2;
	}
}
