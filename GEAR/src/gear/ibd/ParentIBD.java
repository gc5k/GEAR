package gear.ibd;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Map.Entry;

import gear.CmdArgs;
import gear.ConstValues;
import gear.family.pedigree.file.MapFile;
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
	private MapFile snpMap;
	private double[][] ibd;

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
		snpMap = pp.getMapData();
	}

	public void getIBD()
	{
		PrintStream ibdFile = FileUtil.CreatePrintStream(CmdArgs.INSTANCE.out + ".ibd.raw");

		PrintStream ibd2File = FileUtil.CreatePrintStream(CmdArgs.INSTANCE.out + ".ibd");

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

					ibd = new double[per2.getNumMarkers()][2];

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
					
					guessIBD(ibd, 0);
					guessIBD(ibd, 1);
					ibd2File
					.print(per1.getFamilyID() + " " + per1
							.getPersonID() + " " + per2.getFamilyID() + " " + per2
							.getPersonID() + " ");
					for (int k = 0; k < F.getNumMarkers(); k++)
					{
						ibd2File.print(ibd[k][0] + " ");
					}
					ibd2File.println();
			
					ibd2File
					.print(per1.getFamilyID() + " " + per1
					.getPersonID() + " " + per2.getFamilyID() + " " + per2
					.getPersonID() + " ");
					for (int k = 0; k < F.getNumMarkers(); k++)
					{
						ibd2File.print(ibd[k][1] + " ");
					}
					ibd2File.println();

				}
			}
		}
		ibdFile.close();
		ibd2File.close();
	}

	private void guessIBD(double[][]ibd, int pidx)
	{
		int idx0 = 0;
		int idx1 = 0;
		boolean end = false;
		if (Math.abs(ibd[0][pidx] - 0.5) < ConstValues.EPSILON)
		{
			//find the first nonzero
			while (Math.abs(ibd[idx0][pidx] - 0.5) < ConstValues.EPSILON)
			{
				idx0++;
				if (idx0 == ibd.length)
				{
					end = true;
					break;
				}
			}
			if (!end)
			{
				for (int i = 0; i < idx0; i++)
				{
					ibd[i][pidx]= ibd[idx0][pidx];
				}
			}
		}

		//find the first zero
		if (!end)
		{
			while (Math.abs(ibd[idx0][pidx] - 0.5) >= ConstValues.EPSILON)
			{
				idx0++;
				if (idx0 == ibd.length)
				{
					end = true;
					break;
				}
			}
			idx1 = idx0;
			idx0--;
		}

		while (!end)
		{
			//right boundary first nonzero
			while (Math.abs(ibd[idx1][pidx] - 0.5) < ConstValues.EPSILON)
			{
				idx1++;
				if (idx1 == ibd.length)
				{
					end = true;
					break;
				}
			}

			if(end)
			{
				break;
			}

			double pos1 = snpMap.getSNP(idx0).getPosition();
			double pos2 = snpMap.getSNP(idx1).getPosition();
			double w = pos2 - pos1;
			for(int i = idx0+1; i < idx1; i++)
			{
				double pos = snpMap.getSNP(i).getPosition();
				double w1 = pos - pos1;
				double w2 = pos2 - pos;
				ibd[i][pidx] = ibd[idx0][pidx] * w2/w + ibd[idx1][pidx] * w1/w;
			}

			//left boundary first zero
			while(Math.abs(ibd[idx1][pidx] - 0.5) >= ConstValues.EPSILON)
			{
				idx1++;
				if(idx1 == ibd.length)
				{
					end = true;
					break;
				}
			}
			idx0 = idx1;
			idx0--;
		}

		
		if (Math.abs(ibd[ibd.length-1][pidx] - 0.5) < ConstValues.EPSILON)
		{
			boolean head = false;
			idx0 = ibd.length-1;
			while (Math.abs(ibd[idx0][pidx] - 0.5) < ConstValues.EPSILON)
			{
				idx0--;
				if (idx0 == 0)
				{
					head = true;
					break;
				}
			}
			if (!head)
			{
				for(int i = idx0+1; i< ibd.length; i++)
				{
					ibd[i][pidx] = ibd[idx0][pidx];
				}
			}
		}
	}

	/**
	 * 
	 * @param fg father's genotype
	 * @param mg mother's genotype
	 * @param kg1 kid1's genotype
	 * @param kg2 kid2's genotype
	 * @return a two-element array indicating whether the two allele
	 * pairs are IBD, where 0 stands for non-IBD, 0.5 stands for unknown,
	 * and 1 stands for IBD. 
	 */
	private double[] quickIBD(int fg, int mg, int kg1, int kg2) 
	{
		double[] ibd = {0.5, 0.5};
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
					ibd[0] = 0; 
				}
				else if (mg == 1)
				{
					ibd[1] = 0;
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
				ibd[0] = 0;
				ibd[1] = 0;
			}
		}
		return ibd;
	}
}
