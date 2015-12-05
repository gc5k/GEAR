package gear.util.stat;

import gear.util.Logger;

import java.util.ArrayList;
import java.util.Collections;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistributionImpl;

public class PrecisePvalue
{
	private static NormalDistributionImpl unitNormal = new NormalDistributionImpl(0, 1);

	public static double ZcumulativeProbability(double x)
	{
		//http://en.wikipedia.org/wiki/Q-function#cite_note-6
		//IEEE Communications Letters, 2007, 11, 644
		double p1 = (1-Math.exp(-1.4*x))*Math.exp(-x*x/2)/(1.135 * Math.sqrt(2*Math.PI) * x);
		return p1;
	}

	public static double TwoTailZcumulativeProbability(double x)
	{
		//http://en.wikipedia.org/wiki/Q-function#cite_note-6
		//IEEE Communications Letters, 2007, 11, 644
		double p1 = ZcumulativeProbability(x);
		return p1*2;
	}
	
	public static double getPvalue4Z_2Tail(double z)
	{
		double p = 0;
		try
		{
			if (Math.abs(z) < 8)
			{
				p = (1-unitNormal.cumulativeProbability(Math.abs(z)))*2;
			}
			else
			{
				p = PrecisePvalue.TwoTailZcumulativeProbability(Math.abs(z));
			}
		}
		catch (MathException e)
		{
			Logger.printUserError(e.toString());
		}

		return p;

	}
	
	public static double getPvalue4Z_1Tail(double z)
	{
		double p = 0;
		try
		{
			if (Math.abs(z) < 8)
			{
				p = unitNormal.cumulativeProbability(Math.abs(z));
			}
			else
			{
				p = z > 0 ? (1 - PrecisePvalue.ZcumulativeProbability(z)) : PrecisePvalue.ZcumulativeProbability(z);
			}
		}
		catch (MathException e)
		{
			Logger.printUserError(e.toString());
		}

		return p;
	}

}
