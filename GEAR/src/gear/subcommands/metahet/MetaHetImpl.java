package gear.subcommands.metahet;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.BufferedReader;
import gear.util.Logger;

public class MetaHetImpl  extends CommandImpl
{

	@Override
	public void execute(CommandArguments cmdArgs)
	{
		mhArgs = (MetaHetArguments) cmdArgs;
		readXMFile();
		readOMFile();
		if(mhArgs.isQT())
		{
			CalQTHet();
		}
		else
		{
			CalCCHet();
		}
		writeHMat();
	}

	public void CalQTHet()
	{
		hMat = new double[cohortN][cohortN];

		for(int i = 0; i < xMat.length; i++)
		{
			double n_i = oMat[i][i];
			for(int j = i+1; j < xMat[i].length; j++)
			{
				double x = xMat[j][i];
				double kappa = xMat[i][j];
				double n_j = oMat[i][j];
				double on = oMat[j][i];
				double rho = on/Math.sqrt(n_i * n_j);

				double H = x/(mhArgs.getMe() * (1-rho*kappa)) - 1;
				double sd = Math.sqrt(2/Math.sqrt(mhArgs.getMe())) / (1-rho*kappa);
				hMat[i][j] = H/sd;
				hMat[j][i] = H;
			}
		}		
	}
	
	public void CalCCHet()
	{
		double[][] ccSize = mhArgs.getCCsize();
		hMat = new double[cohortN][cohortN];
		if(ccSize.length != cohortN)
		{
			Logger.printUserError("the sample size of the case-control cohort does not match. Please check '" + mhArgs.getCCBatchFile() + ".'");
			System.exit(0);
		}

		for(int i = 0; i < xMat.length; i++)
		{
			double r_i = ccSize[i][0]/ccSize[i][1];
			double n_i = ccSize[i][0]+ccSize[i][1];
			for(int j = i+1; j < xMat[i].length; j++)
			{
				double x = xMat[j][i];
				double kappa = xMat[i][j];
				double r_j = ccSize[j][0]/ccSize[j][1];
				double n_j = ccSize[j][0]+ccSize[j][1];
				double n_cs = oMat[i][j];
				double n_ctrl = oMat[j][i];
				double RR_root = Math.sqrt(r_i * r_j);
				double rho = (n_cs * RR_root + n_ctrl / RR_root)/Math.sqrt(n_i * n_j);
				
				double H = x/(mhArgs.getMe() * (1-rho*kappa)) - 1;
				double sd = Math.sqrt(2/Math.sqrt(mhArgs.getMe())) / (1-rho*kappa);
				hMat[i][j] = H/sd;
				hMat[j][i] = H;
			}
		}
	}

	private void writeHMat()
	{
		PrintWriter hwriter = null;
		try 
		{
			hwriter = new PrintWriter(new BufferedWriter(new FileWriter(mhArgs.getOutRoot() + ".hm")));
		}
		catch (IOException e)
		{
			Logger.handleException(e, "An I/O exception occurred when writing '" + mhArgs.getOutRoot() + ".hm" + "'.");
		}

		for (int i = 0; i < hMat.length; i++)
		{
			for (int j = 0; j < hMat[i].length; j++)
			{
				hwriter.print(String.format("%.4f", hMat[i][j]) + " ");
			}
			hwriter.println();
		}
		hwriter.close();

	}
	
	private void readXMFile()
	{
		BufferedReader bf = BufferedReader.openTextFile(mhArgs.getXMFile(), "cm file.");
		Logger.printUserLog("Reading '" + mhArgs.getXMFile() + "'.");

		String[] d = null;
		int cnt = 0;
		while ( (d = bf.readTokens())!= null )
		{
			if (cnt == 0)
			{
				cohortN = d.length;
				xMat = new double[cohortN][cohortN];
			}

			if (d.length != cohortN )
			{
				Logger.printUserError("incorrect '" + mhArgs.getXMFile() + "'.");
				System.exit(0);
			}

			for (int i = 0; i < d.length; i++)
			{
				xMat[cnt][i] = Double.parseDouble(d[i]);
			}
			cnt++;
		}

		Logger.printUserLog("A " +cohortN + "X" + cohortN + " X statistic matrix is read in.");
	}

	public void readOMFile()
	{
		BufferedReader bf = null;
		if(mhArgs.isQT())
		{
			bf = BufferedReader.openTextFile(mhArgs.getQMFile(), "qm file.");
			Logger.printUserLog("Reading '" + mhArgs.getQMFile() + "'.");
		}
		else
		{
			bf = BufferedReader.openTextFile(mhArgs.getCCMFile(), "ccm file.");
			Logger.printUserLog("Reading '" + mhArgs.getCCMFile() + "'.");			
		}

		String[] d = null;
		int cnt = 0;
		while ( (d = bf.readTokens())!= null )
		{
			if (cnt == 0)
			{
				cohortN = d.length;
				oMat = new double[cohortN][cohortN];
			}

			if (d.length != cohortN )
			{
				String s = mhArgs.isQT() ? mhArgs.getQMFile():mhArgs.getCCMFile();
				Logger.printUserError("incorrect '" + s + "'.");
				System.exit(0);
			}

			for (int i = 0; i < d.length; i++)
			{
				oMat[cnt][i] = Double.parseDouble(d[i]);
			}
			cnt++;
		}
		Logger.printUserLog("A " +cohortN + "X" + cohortN + " overlapping sample matrix is read in.");
	}

	private MetaHetArguments mhArgs;
	private double[][] xMat;
	private double[][] oMat;
	private double[][] hMat;
	private int cohortN;
	
}
