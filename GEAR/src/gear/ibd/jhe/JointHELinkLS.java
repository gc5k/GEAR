package gear.ibd.jhe;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;

import org.apache.commons.math.stat.regression.OLSMultipleLinearRegression;

import gear.CmdArgs;
import gear.data.InputDataSet;
import gear.data.SubjectID;
import gear.util.BufferedReader;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;

public class JointHELinkLS
{
	private String ibdFile;
	private String pheFile;
	private int pheIdx = 0;
	private ArrayList<String> ibdID1 = NewIt.newArrayList();
	private ArrayList<String> ibdID2 = NewIt.newArrayList();
	private double[][] pibd;
	private double[][] mibd;

	private HashMap<String, Integer> pheID = NewIt.newHashMap();
	private HashMap<String, Boolean> keep = NewIt.newHashMap();
	private ArrayList<Integer> keepIBD = NewIt.newArrayList();
	
	private double[][] phe;
	private double[] Y;

	public JointHELinkLS ()
	{
		pheFile = CmdArgs.INSTANCE.getHEArgs().getPheno();
		ibdFile = CmdArgs.INSTANCE.ibdFile;
		pheIdx = CmdArgs.INSTANCE.getHEArgs().getTargetTraitOptionValue() -1;
		initial();
	}

	public void initial()
	{
		readPhenotypes();
		readIBD();
		lineup();
	}

	private void lineup()
	{
		//need to take care of the missing data
		ArrayList<Double> HEphe = NewIt.newArrayList();
		for (int i = 0; i < ibdID1.size(); i++)
		{

			if(keep.containsKey(ibdID1.get(i)) && keep.containsKey(ibdID2.get(i)))
			{
				if (keep.get(ibdID1.get(i)) && keep.get(ibdID2.get(i)))
				{
					int idx1 = pheID.get(ibdID1.get(i)).intValue();
					int idx2 = pheID.get(ibdID2.get(i)).intValue();
					double d= phe[idx1][0] - phe[idx2][0];
					HEphe.add(d*d);
					keepIBD.add(i);
				}
			}
		}

		//make ibd
		double[][] tpibd = new double[keepIBD.size()][pibd[0].length];
		double[][] tmibd = new double[keepIBD.size()][mibd[0].length];
		
		for(int i = 0; i < keepIBD.size(); i++)
		{
			int id = keepIBD.get(i).intValue();
			System.arraycopy(pibd[id], 0, tpibd[id], 0, pibd[id].length);
			System.arraycopy(mibd[id], 0, tmibd[id], 0, mibd[id].length);
		}

		pibd = new double[keepIBD.size()][tpibd[0].length];
		mibd = new double[keepIBD.size()][tmibd[0].length];

		for(int i = 0; i < keepIBD.size(); i++)
		{
			int id = keepIBD.get(i).intValue();
			System.arraycopy(tpibd[id], 0, pibd[id], 0, pibd[id].length);
			System.arraycopy(tmibd[id], 0, mibd[id], 0, mibd[id].length);
		}

		//make Y
		Y = new double[HEphe.size()];

		for (int i = 0; i < Y.length; i++)
		{
			Y[i] = HEphe.get(i);
		}

		Logger.printUserLog(keepIBD.size() + " phenotypes matched with IBD scores.");
	}

	public void JHE()
	{
		JHEpm();
	}

	private void JHEAve()
	{
		
	}

	private void JHEpm()
	{
		StringBuffer sb = new StringBuffer();
		sb.append(CmdArgs.INSTANCE.out);
		sb.append(".helink");
		PrintStream ps = FileUtil.CreatePrintStream(sb.toString());
		OLSMultipleLinearRegression regression = new OLSMultipleLinearRegression();
		Logger.printUserLog("Started scannning...");
		ps.println("mean\tb1(paternal)\tb2(maternal)");
		for (int i = 0; i < pibd[0].length; i++)
		{
			double [][] x = new double[Y.length][2];
			for (int j = 0; j < pibd.length; j++)
			{
				x[j][0] = pibd[j][i];
				x[j][1] = mibd[j][i];
			}
			regression.newSampleData(Y, x);
			double b[] = regression.estimateRegressionParameters();
			
			for (int j = 0; j < b.length; j++)
			{
				ps.print(b[j] + " ");
			}
			ps.println();
		}
		ps.close();
		Logger.printUserLog("Finished HE scanning. The result has been saved in " + sb.toString());
	}
	
	
	private void readPhenotypes()
	{
		InputDataSet ds = new InputDataSet();
		Logger.printUserLog("Reading " + pheFile);
		ds.readPhenotypeFile(pheFile);
		Logger.printUserLog("Read " + ds.getNumberOfSubjects() + " individuals, " + ds.getNumberOfTraits() + " phenotype(s) in " + pheFile);

		if (ds.getNumberOfTraits() < pheIdx)
		{
			Logger.printUserError("the index for the selected phenotype is too large! Only " + ds.getNumberOfTraits() + " phenotypes.");
			System.exit(1);
		}

		phe = new double[ds.getNumberOfSubjects()][1];
		int cn = 0;
		int Cn = 0;
		for (int i = 0; i < ds.getNumberOfSubjects(); i++)
		{
			SubjectID subID = ds.getSubjectID(i);
			StringBuffer sb = new StringBuffer(subID.getFamilyID());
			sb.append(".");
			sb.append(subID.getIndividualID());
			for (int j = 0; j < 1; j++)
			{
				if (ds.isPhenotypeMissing(i, pheIdx))
				{
					keep.put(sb.toString(), false);
					cn++;
				}
				else
				{
					pheID.put(sb.toString(), Cn++);
					keep.put(sb.toString(), true);
					phe[i][j] = ds.getPhenotype(i, j);
				}
			}
		}
		
		Logger.printUserLog(Cn + " nonmissing values in " + pheFile + " for the selected trait.");
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
