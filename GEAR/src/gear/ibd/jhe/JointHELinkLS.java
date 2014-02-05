package gear.ibd.jhe;

import java.util.ArrayList;
import java.util.HashMap;

import gear.CmdArgs;
import gear.data.InputDataSet;
import gear.util.BufferedReader;
import gear.util.Logger;
import gear.util.NewIt;

public class JointHELinkLS
{
	private String ibdFile;
	private String pheFile;
	private ArrayList<String> ibdID1 = NewIt.newArrayList();
	private ArrayList<String> ibdID2 = NewIt.newArrayList();
	private double[][] pibd;
	private double[][] mibd;

	private HashMap<String, Integer> pheID = NewIt.newHashMap();
	private double[][] phe;

	public JointHELinkLS ()
	{
		pheFile = CmdArgs.INSTANCE.getHEArgs().getPheno();
		ibdFile = CmdArgs.INSTANCE.ibdFile;
		initial();
	}
	
	public void initial()
	{
		readPhenotypes();
		readIBD();
	}

	private void readPhenotypes()
	{
		InputDataSet ds = new InputDataSet();
		Logger.printUserLog("Reading " + pheFile);
		ds.readPhenotypeFile(pheFile);

		phe = new double[ds.getNumberOfSubjects()][ds.getNumberOfTraits()];
		for (int i = 0; i < ds.getNumberOfSubjects(); i++)
		{
			StringBuffer sb = new StringBuffer(ds.getSubjectID(i).getFamilyID());
			sb.append(".");
			sb.append(ds.getSubjectID(i).getIndividualID());
			pheID.put(ds.toString(), i);
			for (int j = 0; j < ds.getNumberOfTraits(); j++)
			{
				phe[i][j] = ds.getPhenotype(i, j);
			}
		}
		Logger.printUserLog("Read " + ds.getNumberOfSubjects() + " individuals, " + ds.getNumberOfTraits() + " phenotype(s) in " + pheFile);
	}
	
	private void readIBD()
	{
		BufferedReader reader = BufferedReader.openTextFile(ibdFile, "ibd");

		String[] tokens1 = reader.readTokensAtLeast(4), tokens2;
		
		int numCols = tokens1.length;
		if (numCols % 2 != 0)
		{
			reader.errorPreviousLine("Odd number of columns (" + numCols + " columns).");
		}
		
		ArrayList<double[]> pibd = new ArrayList<double[]>();
		ArrayList<double[]> mibd = new ArrayList<double[]>();

		do
		{	
			double[] ibdThisDad = new double[numCols - 4];
			for (int j = 4; j < numCols; ++j)
			{
				try
				{
					ibdThisDad[j - 4] = Double.parseDouble(tokens1[j]);
				}
				catch (NumberFormatException e)
				{
					reader.errorPreviousLine("'" + tokens1[j] + "' is not a valid floating point number.");
				}
			}
			pibd.add(ibdThisDad);

			tokens2 = reader.readTokens(numCols);
			
			String id11 = tokens1[0] + "." + tokens1[1];
			String id12 = tokens1[2] + "." + tokens1[3];
			String id21 = tokens2[0] + "." + tokens2[1];
			String id22 = tokens2[2] + "." + tokens2[3];

			if (!id11.equals(id21) && !id12.equals(id22))
			{
				reader.errorPreviousLine("The IDs in this line and the above line do not match.");
			}
			
			ibdID1.add(id11);
			ibdID2.add(id12);
			
			double[] ibdThisMom = new double[numCols - 4];
			for (int j = 4; j < numCols; ++j)
			{
				try
				{
					ibdThisMom[j - 4] = Double.parseDouble(tokens2[j]);
				}
				catch (NumberFormatException e)
				{
					reader.errorPreviousLine("'" + tokens2[j] + "' is not a valid floating point number.");
				}
			}
			mibd.add(ibdThisMom);
		} while ((tokens1 = reader.readTokens(numCols)) != null);
		
		this.pibd = pibd.toArray(new double[0][]);
		this.mibd = mibd.toArray(new double[0][]);
		
		Logger.printUserLog("Read " + ibdID1.size() + " pairs in " + ibdFile);
	}
}
