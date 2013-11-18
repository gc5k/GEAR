package gear.metawatchdog.decrypt;

import gear.CmdArgs;
import gear.data.SubjectID;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;

import org.apache.commons.math.stat.regression.*;
import org.apache.commons.math.stat.StatUtils;

public class MetaWatchdog2
{
	private final String delim = "\\s+";

	private double[][] phe1;
	private ArrayList<String> ID1 = NewIt.newArrayList();
	private String f1 = null;
	
	private double[][] phe2;
	private ArrayList<String> ID2 = NewIt.newArrayList();
	private String f2 = null;

	private double cutoff = 0.95;

	public MetaWatchdog2 ()
	{

		f1 = CmdArgs.INSTANCE.set1_file;
		f2 = CmdArgs.INSTANCE.set2_file;
		cutoff = CmdArgs.INSTANCE.dog_cutoff;

		Logger.printUserLog("Set1: " + f1);
		Logger.printUserLog("Set2: " + f2);
		Logger.printUserLog("Cutoff: " + cutoff);
	}

	public void Bark()
	{
		readPhenotypes1(f1);
		readPhenotypes2(f2);

		PrintStream predictorFile = FileUtil.CreatePrintStream(CmdArgs.INSTANCE.out + ".watchdog2");

		int test = phe1[0].length <= phe2[0].length ? phe1[0].length : phe2[0].length; 

		int cnt = 0;
		for (int i = 0; i < phe1.length; i++)
		{
			for (int j = 0; j < phe2.length; j++)
			{
				double[][] dat = new double[test][2];
				for (int k = 0; k < test; k++)
				{
					dat[k][0] = phe1[i][k];
					dat[k][1] = phe2[j][k];
				}
				SimpleRegression sr = new SimpleRegression();
				sr.addData(dat);
				double b = sr.getSlope();

				if(b > cutoff)
				{
					predictorFile.println(ID1.get(i) + "\t" +ID2.get(j) + "\t" + b + "\t" + sr.getSlopeStdErr() + "\t" + sr.getN());
					cnt++;
				}
			}
		}
		predictorFile.close();
		Logger.printUserLog("In total " + phe1.length * phe2.length + " pairs were compared.");
		
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

		Logger.printUserLog("Read " + tmp.size() + " individuals, " + (numCols-2) + " scores in " + fileName);
		Logger.printUserLog("Standardization the scores for each individual in set 1");
		for (int i = 0; i < phe1.length; i++)
		{
			phe1[i] = StatUtils.normalize(phe1[i]);
		}

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

		Logger.printUserLog("Read " + tmp.size() + " individuals, " + (numCols-2) + " scores in " + fileName);
		Logger.printUserLog("Standardization the scores for each individual in set 2");
		for (int i = 0; i < phe2.length; i++)
		{
			phe2[i] = StatUtils.normalize(phe2[i]);
		}				
	}

}
