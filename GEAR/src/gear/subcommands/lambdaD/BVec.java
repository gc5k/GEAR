package gear.subcommands.lambdaD;

import java.util.ArrayList;

import org.apache.commons.math.MathException;
import org.apache.commons.math.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math.stat.regression.SimpleRegression;

import gear.util.Logger;
import gear.util.NewIt;

public class BVec
{
	public BVec()
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
		s.add(b1*b1/(se1*se1));
		s.add(b2*b2/(se2*se2));
		sum.add(s);
	}

	public ArrayList<Double> get(int i)
	{
		return sum.get(i);
	}

	public void setSelected()
	{
		selIdx = new int[sum.size()];
		for(int i = 0; i < selIdx.length; i++)
		{
			selIdx[i] = i;
		}
	}
	
	public void setSelected(int[] idx)
	{
		selIdx = idx;
	}

	public void CalCorrelation()
	{
		double[][] beta = new double[selIdx.length][2];
		double[][] zscore = new double[selIdx.length][2];
		double[][] chiscore = new double[selIdx.length][2];

		double Z12 = 0, Z1sq = 0, Z2sq = 0;

		for(int i = 0; i < selIdx.length; i++)
		{
			ArrayList<Double> s = sum.get(selIdx[i]);
			beta[i][0] = s.get(0);
			beta[i][1] = s.get(1);
			
			zscore[i][0] = s.get(4);
			zscore[i][1] = s.get(5);
			Z12 += zscore[i][0] * zscore[i][1];
			Z1sq += zscore[i][0] * zscore[i][0];
			Z2sq += zscore[i][1] * zscore[i][1];

			chiscore[i][0] = s.get(6);
			chiscore[i][1] = s.get(7);
		}

		PearsonsCorrelation bpc = new PearsonsCorrelation(beta);
		PearsonsCorrelation zpc = new PearsonsCorrelation(zscore);
		PearsonsCorrelation chipc = new PearsonsCorrelation(chiscore);
		
		rb = bpc.getCorrelationMatrix().getData()[1][0];
		rz = zpc.getCorrelationMatrix().getData()[1][0];
		rchi = chipc.getCorrelationMatrix().getData()[1][0];
		
		try
		{
			prb = bpc.getCorrelationPValues().getData()[1][0];
			prz = zpc.getCorrelationPValues().getData()[1][0];
			prchi = chipc.getCorrelationPValues().getData()[1][0];
		}
		catch (MathException e)
		{
			e.printStackTrace();
		}
		
		blm = new SimpleRegression();
		blm.addData(beta);
		
		
		zlm = new SimpleRegression();
		zlm.addData(zscore);		
		
		chilm = new SimpleRegression();
		chilm.addData(chiscore);

		/////////////////
		rg = Z12/Math.sqrt(((Z1sq - zscore.length) * (Z2sq - zscore.length)));
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

	public double getRg()
	{
		return rg;
	}

	public void printOut()
	{
		Logger.printUserLog("===Random effect model===================");
		Logger.printUserLog("Effective number of markers is: " + selIdx.length);
		Logger.printUserLog("Genetic effect correlation: " + rb);
		Logger.printUserLog("p-value for z score (two-tails): " + prb);

		Logger.printUserLog("Effects regression:");
		Logger.printUserLog("Effect regression intercept: " + blm.getIntercept() + " (" + blm.getInterceptStdErr() + ")");
		Logger.printUserLog("Effect regression coefficient: " + blm.getSlope() + " (" + blm.getSlopeStdErr() + ")");

		
		Logger.printUserLog("z score correlation: " + rz);
		Logger.printUserLog("p-value for z score (two-tails): " + prz);

		Logger.printUserLog("Z score regression:");
		Logger.printUserLog("Effect regression intercept: " + zlm.getIntercept() + " (" + zlm.getInterceptStdErr() + ")");
		Logger.printUserLog("Effect regression coefficient: " + zlm.getSlope() + " (" + zlm.getSlopeStdErr() + ")");

		Logger.printUserLog("Chisq score correlation: " + rchi);
		Logger.printUserLog("p-value for z score (two-tails): " + prchi);

		Logger.printUserLog("Chisq score regression:");
		Logger.printUserLog("Effect regression intercept: " + chilm.getIntercept() + " (" + chilm.getInterceptStdErr() + ")");
		Logger.printUserLog("Effect regression coefficient: " + chilm.getSlope() + " (" + chilm.getSlopeStdErr() + ")");
		
		Logger.printUserLog("Genetic correlation(rg):" + rg);
	}

	ArrayList<ArrayList<Double>> sum = null;
	
	private int[] selIdx;
	
	private double rb;
	private double prb;

	private double rz;
	private double prz;

	private double rchi;
	private double prchi;

	private double rg = 0;

	private SimpleRegression blm;
	private SimpleRegression zlm;
	private SimpleRegression chilm;
}
