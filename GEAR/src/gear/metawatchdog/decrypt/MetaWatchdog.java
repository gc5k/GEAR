package gear.metawatchdog.decrypt;

import gear.CmdArgs;
import gear.data.SubjectID;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;
import gear.util.stat.PCA.NCGStatUtils;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;

public class MetaWatchdog
{
	private final String delim = "\\s+";

	private double[][] phe1;
	private ArrayList<String> ID1 = NewIt.newArrayList();
	private String f1 = null;
	
	private double[][] phe2;
	private ArrayList<String> ID2 = NewIt.newArrayList();
	private String f2 = null;

	private double alpha = 0.05;
	private NormalDistribution nd = new NormalDistributionImpl();
	private double T = 0;

	public MetaWatchdog ()
	{

		f1 = CmdArgs.INSTANCE.set1_file;
		f2 = CmdArgs.INSTANCE.set2_file;
		alpha = CmdArgs.INSTANCE.alpha;
		try
		{
			T = nd.inverseCumulativeProbability(alpha/2 + 0.5) * Math.sqrt(2);
		}
		catch (MathException e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		Logger.printUserLog("set1: " + f1);
		Logger.printUserLog("Set2: " + f2);
		Logger.printUserLog("alpha: " + alpha);
		Logger.printUserLog("cutoff: " + T);
	}

	public void Bark()
	{
		readPhenotypes1(f1);
		readPhenotypes2(f2);

		PrintStream predictorFile = FileUtil.CreatePrintStream(CmdArgs.INSTANCE.out + ".watchdog");

		int test = phe1[0].length <= phe2[0].length ? phe1[0].length : phe2[0].length; 

		boolean flag = true;
		
		int cnt = 0;
		for (int i = 0; i < phe1.length; i++)
		{
			for (int j = 0; j < phe2.length; j++)
			{
				flag = true;
				double[] s = new double[test];
				for (int k = 0; k < test; k++)
				{
					s[k] = Math.abs(phe1[i][k] - phe2[j][k]);
					if (s[k] > T)
					{
						flag = false;
					}
				}
				if (flag)
				{
					cnt++;
					predictorFile.print(ID1.get(i) + "\t" +ID2.get(j) + "\t");
					for (int l = 0; l < test; l++)
					{
						if (l < (test - 1) ) 
						{
							predictorFile.print(s[l] + "\t");
						}
						else
						{
							predictorFile.println(s[l]);							
						}
					}
				}
			}
		}
		predictorFile.close();
		Logger.printUserLog("In total " + cnt + " similar pairs were detected.");
	}

	private void readPhenotypes1(String fileName)
	{
		Logger.printUserLog("Reading " + fileName);
		gear.util.BufferedReader reader = gear.util.BufferedReader.openTextFile(fileName, "phenotype");

		HashSet<SubjectID> subjectsRead = new HashSet<SubjectID>();

		String[] tokens = reader.readTokensAtLeast(3);

		if (tokens == null)
		{
			Logger.printUserError("The phenotype file '" + fileName + "' is empty.");
			System.exit(1);
		}

		int numCols = tokens.length;

		ArrayList<String> tmp = NewIt.newArrayList();

		int indIdx = 0;

		while ((tokens = reader.readTokens(numCols)) != null)
		{
			SubjectID subID = new SubjectID(/*famID*/tokens[0], /*indID*/tokens[1]);

			if (subjectsRead.contains(subID))
			{
				Logger.printUserLog("Individual:" + subID.getFamilyID() + " " + subID.getIndividualID() + " at line " + (indIdx+1) + " excluded due to double entry.");
			}
			else
			{
				StringBuffer id = new StringBuffer();
				id.append(subID.getFamilyID() + "\t" + subID.getIndividualID());

				StringBuffer sb = new StringBuffer();
				sb.append(id.toString());
				for (int i = 2; i < numCols; i++)
				{
					sb.append("\t" + tokens[i]);
				}
				tmp.add(sb.toString());
			}
			indIdx++;
		}
		reader.close();

		phe1 = new double[tmp.size()][numCols-2];
		for(int i = 0; i < tmp.size(); i++) 
		{
			String[] tk = tmp.get(i).split(delim);
			StringBuffer id = new StringBuffer();
			id.append(tk[0] + "\t" + tk[1] + "\t");
			ID1.add(id.toString());

			for(int j = 0; j < tk.length - 2; j++)
			{
				phe1[i][j] = Double.parseDouble(tk[j+2]);
			}
		}
		Logger.printUserLog("Read " + tmp.size() + " individuals, " + (numCols-2) + " scores. in " + fileName);
		phe1 = NCGStatUtils.standardize(phe1, true);		
	}

	private void readPhenotypes2(String fileName)
	{
		Logger.printUserLog("Reading " + fileName);
		gear.util.BufferedReader reader = gear.util.BufferedReader.openTextFile(fileName, "phenotype");

		HashSet<SubjectID> subjectsRead = new HashSet<SubjectID>();

		String[] tokens = reader.readTokensAtLeast(3);

		if (tokens == null)
		{
			Logger.printUserError("The phenotype file '" + fileName + "' is empty.");
			System.exit(1);
		}

		int numCols = tokens.length;

		ArrayList<String> tmp = NewIt.newArrayList();

		int indIdx = 0;

		while ((tokens = reader.readTokens(numCols)) != null)
		{
			SubjectID subID = new SubjectID(/*famID*/tokens[0], /*indID*/tokens[1]);

			if (subjectsRead.contains(subID))
			{
				Logger.printUserLog("Individual:" + subID.getFamilyID() + " " + subID.getIndividualID() + " at line " + (indIdx+1) + " excluded due to double entry.");
			}
			else
			{
				StringBuffer id = new StringBuffer();
				id.append(subID.getFamilyID() + "\t" + subID.getIndividualID());

				StringBuffer sb = new StringBuffer();
				sb.append(id.toString());
				for (int i = 2; i < numCols; i++)
				{
					sb.append("\t" + tokens[i]);
				}
				tmp.add(sb.toString());
			}
			indIdx++;
		}

		reader.close();

		phe2 = new double[tmp.size()][numCols-2];
		for(int i = 0; i < tmp.size(); i++) 
		{
			String[] tk = tmp.get(i).split(delim);
			StringBuffer id = new StringBuffer();
			id.append(tk[0] + "\t" + tk[1] + "\t");
			ID2.add(id.toString());

			for(int j = 0; j < tk.length - 2; j++)
			{
				phe2[i][j] = Double.parseDouble(tk[j+2]);
			}
		}

		Logger.printUserLog("Read " + tmp.size() + " individuals, " + (numCols-2) + " scores. in " + fileName);
		phe2 = NCGStatUtils.standardize(phe2, true);
	}

}
