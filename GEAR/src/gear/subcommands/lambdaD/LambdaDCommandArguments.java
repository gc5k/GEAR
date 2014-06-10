package gear.subcommands.lambdaD;

import gear.subcommands.CommandArguments;
import gear.util.FileUtil;
import gear.util.Logger;

public class LambdaDCommandArguments extends CommandArguments
{

	public void setMeta1(String m1)
	{
		FileUtil.exists(m1);
		this.m1 = m1;
	}

	public String getMeta1()
	{
		return m1;
	}
	
	public void setMeta2(String m2)
	{
		FileUtil.exists(m2);
		this.m2 = m2;
	}

	public String getMeta2()
	{
		return m2;
	}

	public void setCC(String[] cc)
	{
		ccSize[0] = Double.parseDouble(cc[0]);
		ccSize[1] = Double.parseDouble(cc[1]);
		ccSize[2] = Double.parseDouble(cc[2]);
		ccSize[3] = Double.parseDouble(cc[3]);

		for (int i = 0; i < ccSize.length; i++)
		{
			if(ccSize[i] <= 1)
			{
				Logger.printUserError("The sample size should be greater than 1.");
				System.exit(0);
			}
		}
		isQT = false;
	}
	
	public double[] getCCsize()
	{
		return ccSize;
	}

	public void setQT(String[] qt)
	{
		qtSize[0] = Double.parseDouble(qt[0]);
		qtSize[1] = Double.parseDouble(qt[1]);
		
		for (int i = 0; i < qtSize.length; i++)
		{
			if (qtSize[i] <= 1)
			{
				Logger.printUserError("The sample size should be greater than 1.");
				System.exit(0);
			}
		}
		isQT = true;
	}

	public double[] getQTsize()
	{
		return qtSize;
	}

	public boolean isQT() 
	{
		return isQT;
	}

	private String m1;
	private String m2;

	private boolean isQT = true;
	private double[] qtSize = {0,0};
	private double[] ccSize = {0,0,0,0};
}
