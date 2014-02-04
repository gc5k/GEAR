package gear.ibd.ibd;

import java.io.PrintStream;

import gear.CmdArgs;
import gear.ConstValues;
import gear.data.Person;
import gear.data.Family;
import gear.data.UniqueRecordSet;
import gear.family.pedigree.file.MapFile;
import gear.family.pedigree.file.PedigreeFile;
import gear.family.plink.PLINKBinaryParser;
import gear.family.plink.PLINKParser;
import gear.util.FileUtil;
import gear.util.Logger;

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
		Family family;

		UniqueRecordSet<Family> families = PedData.getFamilies();
		for (int famIdx = 0; famIdx < families.size(); ++famIdx)
		{
			family = families.get(famIdx);

			for (int personIdx1 = 0; personIdx1 < family.size() - 1; ++personIdx1)
			{
				Person person1 = family.getPerson(personIdx1);
				Person dad = family.getPerson(person1.getDadID());
				Person mom = family.getPerson(person1.getMomID());
				
				if (dad == null || mom == null)
				{
					continue;
				}

				for (int personIdx2 = personIdx1 + 1; personIdx2 < family.size(); ++personIdx2)
				{
					Person person2 = family.getPerson(personIdx2);

					if (!dad.getID().equals(person2.getDadID()) ||
						!mom.getID().equals(person2.getMomID()))
					{
						continue;
					}

					ibd = new double[person2.getNumMarkers()][2];
					
					for (int k = 0; k < dad.getNumMarkers(); k++)
					{
						ibd[k] = quickIBD(dad.getGenotypeScore(k),
										  mom.getGenotypeScore(k),
										  person1.getGenotypeScore(k),
										  person2.getGenotypeScore(k));
					}
					
					
					ibdFile.print(person1.getFamilyID() + " " + person1
									.getPersonID() + " " + person2.getFamilyID() + " " + person2
									.getPersonID() + " ");
					for (int k = 0; k < dad.getNumMarkers(); k++)
					{
						ibdFile.print(ibd[k][0] + " ");
					}
					ibdFile.println();
					
					ibdFile
					.print(person1.getFamilyID() + " " + person1
							.getPersonID() + " " + person2.getFamilyID() + " " + person2
							.getPersonID() + " ");
					for (int k = 0; k < dad.getNumMarkers(); k++)
					{
						ibdFile.print(ibd[k][1] + " ");
					}
					ibdFile.println();
					
					guessIBD(ibd, 0);
					guessIBD(ibd, 1);
					ibd2File
					.print(person1.getFamilyID() + " " + person1
							.getPersonID() + " " + person2.getFamilyID() + " " + person2
							.getPersonID() + " ");
					for (int k = 0; k < dad.getNumMarkers(); k++)
					{
						ibd2File.print(ibd[k][0] + " ");
					}
					ibd2File.println();
			
					ibd2File
					.print(person1.getFamilyID() + " " + person1
					.getPersonID() + " " + person2.getFamilyID() + " " + person2
					.getPersonID() + " ");
					for (int k = 0; k < dad.getNumMarkers(); k++)
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
