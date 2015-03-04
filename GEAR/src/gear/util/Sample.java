package gear.util;

import gear.CmdArgs;

import java.util.Random;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ChiSquaredDistributionImpl;

/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */

public class Sample
{
	public static Random U = new Random(CmdArgs.INSTANCE.simuSeed);

	public static void setSeed(long s)
	{
		U.setSeed(s);
	}

	public static int[] sample(int[] s)
	{
		int[] seq = new int[s.length];
		System.arraycopy(s, 0, seq, 0, s.length);
		int[] index = new int[s.length];
		int len = s.length;
		for (int i = 0; i < index.length; i++)
		{
			int j = U.nextInt(len);
			len--;
			index[i] = seq[j];
			seq[j] = seq[len];
			seq[len] = index[i];
		}
		return index;
	}

	public static int[] SampleIndex(int start, int end, int n)
	{
		int[] seq = new int[end - start + 1];
		int[] index = new int[n];
		for (int i = 0; i < seq.length; i++)
			seq[i] = i + start;
		int len = seq.length;
		for (int i = 0; i < index.length; i++)
		{
			int j = U.nextInt(len);
			len--;
			index[i] = seq[j];
			seq[j] = seq[len];
			seq[len] = index[i];
		}
		return index;
	}

	public static int[] SampleIndexWithReplacement(int start, int end, int n)
	{
		int[] seq = new int[end - start + 1];
		int[] index = new int[n];
		for (int i = 0; i < seq.length; i++)
			seq[i] = i + start;
		int len = seq.length;
		for (int i = 0; i < index.length; i++)
		{
			int j = U.nextInt(len);
			index[i] = seq[j];
		}
		return index;
	}
	
	public static double[] sampleChisq(int Len, int df)
	{
		ChiSquaredDistributionImpl chiDis = new ChiSquaredDistributionImpl(df);
		double[] ChiExp = new double[Len];
		for (int i = 0; i < Len; i++)
		{
			try
			{
				ChiExp[i] = chiDis
						.inverseCumulativeProbability((i + 1) / (Len + 0.05));
				;
			}
			catch (MathException e)
			{
				e.printStackTrace();
			}
		}
		return ChiExp;
	}


}
