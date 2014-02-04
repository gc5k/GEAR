package gear.ibd.jhe;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import gear.CmdArgs;
import gear.data.InputDataSet;
import gear.util.FileUtil;
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
		BufferedReader bf = FileUtil.FileOpen(ibdFile);
		String line = null;
		ArrayList<String> tibd = NewIt.newArrayList();

		int cn = 0;
		int Len = 0;
		boolean quit = false;
		try
		{
			while((line = bf.readLine()) !=null)
			{
				String[] s = line.split("\\s+");
				if (cn == 0)
				{
					Len = s.length;
				}
				else
				{
					if(s.length != Len)
					{
						Logger.printUserError("Line " + (cn+1) + " does not have same number of elements as the first line");
						quit = true;
					}
				}
				cn++;
				tibd.add(line);
			}
		}
		catch (IOException e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		if(quit)
		{
			System.exit(1);
		}
		
		if(tibd.size() % 2 != 0)
		{
			Logger.printUserLog("Odd number of lines in " + ibdFile + ".");
			Logger.printUserLog("gear quitted");
		}

		pibd = new double[tibd.size()][Len-4];
		mibd = new double[tibd.size()][Len-4];

		for (int i = 0; i < tibd.size() /2 ; i++)
		{
			String s1[] = tibd.get(i*2).split("\\s+");
			String s2[] = tibd.get(i*2+1).split("\\s+");
			StringBuffer sb1 = new StringBuffer();
			sb1.append(s1[0]);
			sb1.append(".");
			sb1.append(s1[1]);
			ibdID1.add(sb1.toString());
			
			StringBuffer sb2 = new StringBuffer();
			sb2.append(s1[2]);
			sb2.append(".");
			sb2.append(s1[3]);
			ibdID2.add(sb2.toString());

			StringBuffer sb3 = new StringBuffer();
			sb3.append(s2[0]);
			sb3.append(".");
			sb3.append(s2[1]);
			
			StringBuffer sb4 = new StringBuffer();
			sb4.append(s2[2]);
			sb4.append(".");
			sb4.append(s2[3]);
			if (ibdID1.get(i).compareTo(sb3.toString())!= 0 && ibdID2.get(i).compareTo(sb4.toString()) != 0)
			{
				Logger.printUserError("Line " + (i*2) + " and line " + (i*2 +1) + " do not match ids.");
				Logger.printUserLog("gear quitted.");
				System.exit(1);
			}

			for (int j = 4; j < s1.length; j++)
			{
				pibd[i][j-4] = Double.parseDouble(s1[j]);
			}
			for (int j = 4; j < s2.length; j++)
			{
				mibd[i][j-4] = Double.parseDouble(s2[j]);
			}

		}
		Logger.printUserLog("Read " + ibdID1.size() + " pairs in " + ibdFile);
	}
}
