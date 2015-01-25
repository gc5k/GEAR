package gear.subcommands.randmeta;

import java.util.ArrayList;

import org.apache.commons.math.MathException;
import org.apache.commons.math.stat.correlation.PearsonsCorrelation;

import gear.util.Logger;
import gear.util.NewIt;

public class BetaVec
{
	public BetaVec()
	{
		sum= NewIt.newArrayList();
	}

	public void addStats(double b1, double b2, double se1, double se2)
	{
		ArrayList<Double> s = NewIt.newArrayList();
		s.add(b1);
		s.add(b2);
		s.add(se1);
		s.add(se2);
		s.add(b1/se1);
		s.add(b2/se2);
		sum.add(s);
	}

	public void setSelected()
	{
		selIdx = new int[sum.size()];
	}
	
	public void setSelected(int[] idx)
	{
		selIdx = idx;
	}

	public void CalCorrelation()
	{
		double[][] beta = new double[selIdx.length][2];
		double[][] zscore = new double[selIdx.length][2];
		
		for(int i = 0; i < selIdx.length; i++)
		{
			ArrayList<Double> s = sum.get(selIdx[i]);
			beta[i][0] = s.get(0);
			beta[i][1] = s.get(1);
			
			zscore[i][0] = s.get(4);
			zscore[i][1] = s.get(5);
		}

		PearsonsCorrelation bpc = new PearsonsCorrelation(beta);
		PearsonsCorrelation zpc = new PearsonsCorrelation(zscore);
		
		rb = bpc.getCorrelationMatrix().getData()[1][0];
		rz = zpc.getCorrelationMatrix().getData()[1][0];
		
		try
		{
			prb = bpc.getCorrelationPValues().getData()[1][0];
			prz = zpc.getCorrelationPValues().getData()[1][0];
		}
		catch (MathException e)
		{
			e.printStackTrace();
		}
		
	}

	public double getBcorrelation()
	{
		return rb;
	}
	
	public double getZcorrelation()
	{
		return rz;
	}
	
	public double getPvBcorrelation()
	{
		return prb;
	}
	
	public double getPvZcorrelation()
	{
		return prz;
	}

	public void printOut()
	{
		Logger.printUserLog("Effective number of markers is: " + selIdx.length);
		Logger.printUserLog("Genetic effect correlation: " + rb);
		Logger.printUserLog("p-value for z score (two-tails): " + prb);
		Logger.printUserLog("z score correlation: " + rz);
		Logger.printUserLog("p-value for z score (two-tails): " + prz);
	}

	ArrayList<ArrayList<Double>> sum = null;
	
	private int[] selIdx;
	
	private double rb;
	private double prb;

	private double rz;
	private double prz;

}
