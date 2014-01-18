package gear.ibd;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Map.Entry;

import gear.CmdArgs;
import gear.family.pedigree.file.PedigreeFile;
import gear.family.pedigree.genotype.BFamilyStruct;
import gear.family.pedigree.genotype.BPerson;
import gear.family.plink.PLINKBinaryParser;
import gear.family.plink.PLINKParser;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;

public class ParentIBD
{
	private PedigreeFile PedData;
	
	private int[][] ibd;

	public ParentIBD() 
	{
		PLINKParser pp = null;
		if (CmdArgs.INSTANCE.getFileArgs().isSet())
		{
			pp = new PLINKParser(CmdArgs.INSTANCE.getFileArgs()
					.getPed(), CmdArgs.INSTANCE.getFileArgs()
					.getMap());
		}
		else if (CmdArgs.INSTANCE.getBFileArgs(0).isSet())
		{
			pp = new PLINKBinaryParser(CmdArgs.INSTANCE.getBFileArgs(0)
					.getBed(), CmdArgs.INSTANCE.getBFileArgs(0)
					.getBim(), CmdArgs.INSTANCE.getBFileArgs(0)
					.getFam());
		}
		else
		{
			Logger.printUserError("No input files.");
			System.exit(1);
		}
		pp.Parse();
		PedData = pp.getPedigreeData();

	}

	public void getIBD()
	{
		PrintStream ibdFile = FileUtil.CreatePrintStream(CmdArgs.INSTANCE.out + ".ibd");
		
		Enumeration<String> perList1;
		BFamilyStruct fam;

		Hashtable<String, BFamilyStruct> Fam = PedData.getFamilyStruct();
		for (Entry<String, BFamilyStruct> entry : Fam.entrySet())
		{
			fam = (BFamilyStruct) entry.getValue();
			perList1 = fam.getPersonList();

			ArrayList<String> perID = NewIt.newArrayList();
			while (perList1.hasMoreElements())
			{
				String iid1 = perList1.nextElement();
				if (!fam.hasAncestor(iid1))
				{
					continue;
				}
				perID.add(iid1);
			}

			for (int i = 0; i < perID.size(); i++)
			{
				String iid1 = perID.get(i);
				BPerson per1 = fam.getPerson(iid1);

				String fid1 = per1.getDadID();
				String mid1 = per1.getMomID();

				for (int j = i + 1; j < perID.size(); j++)
				{
					String iid2 = perID.get(j);
					BPerson per2 = fam.getPerson(iid2);

					String fid2 = per2.getDadID();
					String mid2 = per2.getMomID();

					if (fid1.compareTo(fid2) != 0 || mid1.compareTo(mid2) != 0)
					{
						continue;
					}

					ibd = new int[per2.getNumMarkers()][2];

					BPerson F = fam.getPerson(fid2);
					BPerson M = fam.getPerson(mid2);

					for (int k = 0; k < F.getNumMarkers(); k++)
					{
						ibd[k] = quickIBD(	F.getGenotypeScore(k),
											M.getGenotypeScore(k),
											per1.getGenotypeScore(k),
											per2.getGenotypeScore(k));
					}

					ibdFile
							.print(per1.getFamilyID() + " " + per1
									.getPersonID() + " " + per2.getFamilyID() + " " + per2
									.getPersonID() + " ");
					for (int k = 0; k < F.getNumMarkers(); k++)
					{
						ibdFile.print(ibd[k][0] + " ");
					}
					ibdFile.println();
					
					ibdFile
					.print(per1.getFamilyID() + " " + per1
							.getPersonID() + " " + per2.getFamilyID() + " " + per2
							.getPersonID() + " ");
					for (int k = 0; k < F.getNumMarkers(); k++)
					{
						ibdFile.print(ibd[k][1] + " ");
					}
					ibdFile.println();
				}
			}
		}
		ibdFile.close();
	}

	private int[] quickIBD(int fg, int mg, int kg1, int kg2) 
	{
		int[] ibd = {0, 0};
		if (fg == 3 || mg == 3 || kg1 == 3 || kg2 ==3)
		{
			return ibd;
		}

		if (Math.abs(fg - mg) == 1)
		{
			if (kg1 == kg2)
			{
				if (kg1 == 1 && kg2 == 1)
				{
					if (fg == 1)
					{
						ibd[0] = 1;
					}
					else if (mg == 1)
					{
						ibd[1] = 1;
					}
				}
				else if (kg1 != 1 && kg2 != 1)
				{
					if (fg == 1)
					{
						ibd[0] = 1;
					}
					else if (mg == 1)
					{
						ibd[1] = 1;
					}
				}
			}
			else
			{
				if (fg == 1)
				{
					ibd[0] = -1; 
				}
				else if (mg == 1)
				{
					ibd[1] = -1;
				}
			}
		}

		if(fg == 1 && mg == 1)
		{
			if (kg1 == 0 && kg2 == 0)
			{
				ibd[0] = 1;
				ibd[1] = 1;
			}
			else if (kg1 == 2 && kg2 == 2)
			{
				ibd[0] = 1;
				ibd[1] = 1;
			}
			else if ( Math.abs(kg1 - kg2 ) == 2)
			{
				ibd[0] = -1;
				ibd[1] = -1;
			}
		}
		return ibd;
	}
}
