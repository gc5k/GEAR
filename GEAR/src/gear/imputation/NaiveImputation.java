package gear.imputation;

import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import gear.CmdArgs;
import gear.ConstValues;
import gear.data.Person;
import gear.family.pedigree.PersonIndex;
import gear.family.pedigree.file.SNP;
import gear.family.plink.PLINKBinaryParser;
import gear.family.plink.PLINKParser;
import gear.family.popstat.GenotypeMatrix;
import gear.family.qc.rowqc.SampleFilter;
import gear.util.Logger;
import gear.util.pop.PopStat;

public class NaiveImputation
{
	private SampleFilter sf;
	private GenotypeMatrix pGM;

	public NaiveImputation()
	{
		initial();
	}

	private void initial ()
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

		sf = new SampleFilter(pp.getPedigreeData());
		pGM = new GenotypeMatrix(sf.getSample(), pp.getMapData());
	}

	public void Imputation() 
	{
		PopStat.Imputation(pGM);
		if (CmdArgs.INSTANCE.makebedFlag)
		{
			writeBFile();
		}
	}

	public void writeBFile()
	{
		DataOutputStream bedout = null;
		PrintWriter fam = null;
		PrintWriter bim = null;
		try
		{
			bedout = new DataOutputStream(new FileOutputStream(CmdArgs.INSTANCE.out + ".bed"));

			fam = new PrintWriter(new BufferedWriter(new FileWriter(CmdArgs.INSTANCE.out
					+ ".fam")));
			bim = new PrintWriter(new BufferedWriter(new FileWriter(CmdArgs.INSTANCE.out
					+ ".bim")));
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}

		ArrayList<PersonIndex> pi = pGM.getSample();

		for (int i = 0; i < pi.size(); i++)
		{
			PersonIndex p = pi.get(i);
			Person bp = p.getPerson();
			fam.println(bp.getFamilyID() + "\t" + bp.getPersonID() + "\t" + bp.getDadID() + "\t" + bp.getMomID() + "\t" + bp.getGender() + "\t" + bp.getAffectedStatus());
		}
		fam.close();

		ArrayList<SNP> snpList = pGM.getSNPList();
		for (int i = 0; i < snpList.size(); i++)
		{
			SNP snp = snpList.get(i);
			bim.println(snp.getChromosome() + "\t" + snp.getName() + "\t" + snp.getDistance() + "\t" + snp.getPosition() + "\t" + snp.getFirstAllele() + "\t" + snp.getSecAllele());
		}
		bim.close();

		try
		{
			bedout.writeByte(ConstValues.PLINK_BED_BYTE1);
			bedout.writeByte(ConstValues.PLINK_BED_BYTE2);
			bedout.writeByte(ConstValues.PLINK_BED_BYTE3);

			for (int i = 0; i < pGM.getNumMarker(); i++)
			{
				byte gbyte = 0;
				int idx = 0;

				for (int j = 0; j < pGM.getNumIndivdial(); j++)
				{

					byte g = pGM.getOriginalGenotypeScore(j, i);
					g <<= 2 * idx;
					gbyte |= g;
					idx++;

					if (j != (pGM.getNumIndivdial() - 1))
					{
						if (idx == 4)
						{
							bedout.writeByte(gbyte);
							gbyte = 0;
							idx = 0;
						}
					}
					else
					{
						bedout.writeByte(gbyte);
					}
				}
			}
			bedout.close();
		}
		catch (IOException e)
		{
			Logger.handleException(e, "An exception occurred when wrtinging bed file.");
		}
	}
}
