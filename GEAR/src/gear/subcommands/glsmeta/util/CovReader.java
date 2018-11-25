package gear.subcommands.glsmeta.util;

import java.util.ArrayList;

import gear.util.BufferedReader;
import gear.util.Logger;
import gear.util.NewIt;

public class CovReader
{
	public CovReader(String covar, int[] covar_idx, int NumCohort)
	{
		setCovFile(covar, NumCohort);

		if (covar != null)
		{
			if (covar_idx == null)
			{
				covar_idx = new int[nCol];
				for(int i = 0; i < nCol; i++)
				{
					covar_idx[i] = i;
				}
			}
			else
			{
				for (int i = 0; i < covar_idx.length; i++)
				{
					if ( (covar_idx[i]-1) >= nCol)
					{
						Logger.printUserLog(covar_idx[i] + " is a too big index for '" + covar + "'.");
						Logger.printUserLog("GEAR quit.");
						System.exit(0);
					}
				}
			}
			covTable = new double[cov.size()][covar_idx.length];
			for(int i = 0; i < covTable.length; i++)
			{
				for(int j = 0; j < covTable[i].length; j++)
				{
					covTable[i][j] = Double.parseDouble(cov.get(i).get(covar_idx[j]-1));
				}
			}
		}


	}

	private void setCovFile(String covar, int NumCohort)
	{
		if (covar != null)
		{
			this.covar = covar;
			BufferedReader reader = BufferedReader.openTextFile(this.covar, "Covar for covariates");
			cov = NewIt.newArrayList();

			String[] tokens = null;
//			int cnt = 0;
			tokens = reader.readTokens();
			ArrayList<String> s = NewIt.newArrayList();
			for(int i = 0; i < tokens.length; i++) s.add(tokens[i]);
			cov.add(s);
			nCol = s.size();

			Logger.printUserLog("'" + covar + "' has " + nCol + " columns.");
			while( (tokens=reader.readTokensAtLeast(nCol)) != null)
			{
				s = NewIt.newArrayList();
				for(int i = 0; i < tokens.length; i++) s.add(tokens[i]);
				cov.add(s);
			}

			if (cov.size() != NumCohort)
			{
				Logger.printUserLog("The number of lines do not match cohort size.");
				Logger.printUserLog("GEAR quit.");
				System.exit(0);
			}
			Logger.printUserLog("Read " + cov.size() + " lines in '" + this.covar + ".'");
		}
	}
	
	public double[][] getCovTable()
	{
		return covTable;
	}

	private String covar = null;
	private ArrayList<ArrayList<String>> cov;
//	private int[] covar_idx;

	private int nCol = 0;
	private double[][] covTable;
}
