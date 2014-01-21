package gear.write;

import java.io.DataOutputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;

import gear.CmdArgs;
import gear.ConstValues;
import gear.data.SubjectID;
import gear.family.pedigree.PersonIndex;
import gear.family.pedigree.file.SNP;
import gear.family.pedigree.genotype.BPerson;
import gear.family.plink.PLINKBinaryParser;
import gear.family.plink.PLINKParser;
import gear.family.qc.rowqc.SampleFilter;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;

public class WriteBedSNPMajor
{

	private ArrayList<PersonIndex> PersonTable;
	private ArrayList<PersonIndex> PrintPerson = NewIt.newArrayList();
	private DataOutputStream os = null;
	private ArrayList<SNP> snpList;

	public WriteBedSNPMajor()
	{
		PLINKParser pp = null;
		if (CmdArgs.INSTANCE.getFileArgs().isSet())
		{
			pp = new PLINKParser(CmdArgs.INSTANCE.getFileArgs()
					.getPed(), CmdArgs.INSTANCE.getFileArgs()
					.getMap());
		} else if (CmdArgs.INSTANCE.getBFileArgs(0).isSet())
		{
			pp = new PLINKBinaryParser(CmdArgs.INSTANCE.getBFileArgs(0)
					.getBed(), CmdArgs.INSTANCE.getBFileArgs(0)
					.getBim(), CmdArgs.INSTANCE.getBFileArgs(0)
					.getFam());
		} else
		{
			Logger.printUserError("No input files.");
			System.exit(1);
		}
		pp.Parse();
		SampleFilter sf = new SampleFilter(pp.getPedigreeData(),
				pp.getMapData());

		snpList = pp.getMapData().getMarkerList();
		PersonTable = sf.getSample();
	}

	public WriteBedSNPMajor(ArrayList<PersonIndex> pt, ArrayList<SNP> sl)
	{
		snpList = sl;
		PersonTable = pt;
	}

	public void WriteFile()
	{
		if(CmdArgs.INSTANCE.orderindFlag)
		{
			OrderInd();
		}
		else
		{
			PrintPerson = PersonTable;
		}
		WriteFile(CmdArgs.INSTANCE.out);
	}

	private void OrderInd()
	{
		Logger.printUserLog("Reading " + CmdArgs.INSTANCE.orderindFile);
		gear.util.BufferedReader reader = gear.util.BufferedReader.openTextFile(CmdArgs.INSTANCE.orderindFile, "phenotype");

		HashSet<SubjectID> subjectsRead = new HashSet<SubjectID>();

		ArrayList<SubjectID> subjectIDs = NewIt.newArrayList();

		String[] tokens = null; //reader.readTokensAtLeast(2);
		int minNumCols = 2;

		int indIdx = 0;

		while ((tokens = reader.readTokens()) != null)
		{
			if (tokens.length < minNumCols)
			{
				Logger.printUserError("There should be at least " + minNumCols + " columns.");
			}

			SubjectID subID = new SubjectID(/*famID*/tokens[0], /*indID*/tokens[1]);

			if (subjectsRead.contains(subID))
			{
				Logger.printUserLog("Individual:" + subID.getFamilyID() + " " + subID.getIndividualID() + " at line " + (indIdx+1) + " excluded due to double entry.");
			}
			else
			{
				StringBuffer id = new StringBuffer();
				id.append(subID.getFamilyID() + "\t" + subID.getIndividualID());

				subjectIDs.add(subID);
				subjectsRead.add(subID);
			}
			indIdx++;
		}
		reader.close();

		Logger.printUserLog("Read " + subjectIDs.size() + " unique individuals in " + CmdArgs.INSTANCE.orderindFile);

		int cn = 0;
		int[] index = new int[subjectIDs.size()];
		Arrays.fill(index, -1);
		for (int i = 0; i < subjectIDs.size(); i++)
		{
			SubjectID subID = subjectIDs.get(i);
			for (int j = 0; j < PersonTable.size(); j++)
			{
				PersonIndex pi = PersonTable.get(j);
				if (subID.getFamilyID().compareTo(pi.getFamilyID()) == 0 && subID.getIndividualID().compareTo(pi.getIndividualID()) == 0)
				{
					cn++;
					index[i] = j;
					break;
				}
			}
		}

		if(cn == 0)
		{
			Logger.printUserLog("no individuals were matched between fam and " + CmdArgs.INSTANCE.orderindFile + ".");
			System.exit(1);
		}

		for(int i = 0; i < index.length; i++)
		{
			if (index[i]!= -1)
			{
				PrintPerson.add(PersonTable.get(index[i]));
			}
		}
		Logger.printUserLog("Matched " + PrintPerson.size() + " individuals in fam");		
	}

	public void WriteFile(String out)
	{
		StringBuffer sbim = new StringBuffer();
		sbim.append(out);
		sbim.append(".bim");
		PrintStream pbim = FileUtil.CreatePrintStream(sbim.toString());
		for (Iterator<SNP> e = snpList.iterator(); e.hasNext();)
		{
			SNP snp = e.next();
			pbim.append(snp.getChromosome() + "\t" + snp.getName() + "\t"
					+ snp.getDistance() + "\t" + snp.getPosition() + "\t"
					+ snp.getFirstAllele() + "\t" + snp.getSecAllele() + "\n");
		}
		pbim.close();

		StringBuffer sfam = new StringBuffer();
		sfam.append(out);
		sfam.append(".fam");
		PrintStream pfam = FileUtil.CreatePrintStream(sfam.toString());
		
		for (int i = 0; i < PrintPerson.size(); i++)
		{
			PersonIndex per = PrintPerson.get(i);
			BPerson bp = per.getPerson();
			pfam.append(bp.getFamilyID() + "\t" + bp.getPersonID() + "\t"
					+ bp.getDadID() + "\t" + bp.getMomID() + "\t"
					+ bp.getGender() + "\t" + bp.getAffectedStatus() + "\n");
		}
		pfam.close();

		StringBuffer sbed = new StringBuffer();
		sbed.append(out);
		sbed.append(".bed");
		try
		{
			os = new DataOutputStream(new FileOutputStream(sbed.toString()));
		} catch (FileNotFoundException e)
		{
			Logger.printUserError("Cannot create file '" + sbed.toString()
					+ "'.");
			Logger.printUserError("Exception Message: " + e.getMessage());
			System.exit(1);
		}

		try
		{
			os.writeByte(ConstValues.PLINK_BED_BYTE1);
			os.writeByte(ConstValues.PLINK_BED_BYTE2);
			os.writeByte(ConstValues.PLINK_BED_BYTE3);

			for (int i = 0; i < snpList.size(); i++)
			{
				byte gbyte = 0;
				int idx = 0;

				int posByte = i >> BPerson.shift;
				int posBit = (i & 0xf) << 1;


				for (int j = 0; j < PrintPerson.size(); j++)
				{
					PersonIndex pi = PrintPerson.get(j);
					BPerson bp = pi.getPerson();
					byte g = bp.getOriginalGenotypeScore(posByte, posBit);

					g <<= 2 * idx;
					gbyte |= g;
					idx++;

					if (j != (PrintPerson.size() - 1))
					{
						if (idx == 4)
						{
							os.writeByte(gbyte);
							gbyte = 0;
							idx = 0;
						}
					}
					else
					{
						os.writeByte(gbyte);
					}
				}
			}
			os.close();
		} catch (IOException e)
		{
			Logger.handleException(
					e,
					"An exception occurred when writing the bed file '"
							+ sbed.toString() + "'.");
		}
	}

}
