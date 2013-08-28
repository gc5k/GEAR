package gear.imputation;

import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import gear.CmdArgs;
import gear.family.pedigree.PersonIndex;
import gear.family.pedigree.file.MapFile;
import gear.family.pedigree.file.SNP;
import gear.family.pedigree.genotype.BPerson;
import gear.family.plink.PLINKBinaryParser;
import gear.family.plink.PLINKParser;
import gear.family.popstat.GenotypeMatrix;
import gear.family.qc.rowqc.SampleFilter;
import gear.sumstat.qc.rowqc.SumStatQC;
import gear.util.Logger;
import gear.util.pop.PopStat;

public class NaiveImputation
{
	private byte byte1 = 108;
	private byte byte2 = 27;
	private byte byte3 = 1;

	private MapFile mapFile;
	private SampleFilter sf;
	private SumStatQC ssQC;
	GenotypeMatrix gm;

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

		sf = new SampleFilter(pp.getPedigreeData(),
				pp.getMapData());
		ssQC = new SumStatQC(pp.getPedigreeData(), pp.getMapData(),
				sf);

		mapFile = ssQC.getMapFile();
		gm = new GenotypeMatrix(ssQC.getSample());
	}

	public void Imputation() 
	{
		PopStat.Imputation(gm);
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

		ArrayList<PersonIndex> pi = ssQC.getSample();

		for (int i = 0; i < pi.size(); i++)
		{
			PersonIndex p = pi.get(i);
			BPerson bp = p.getPerson();
			fam.println(bp.getFamilyID() + "\t" + bp.getPersonID() + "\t" + bp.getDadID() + "\t" + bp.getMomID() + "\t" + bp.getGender() + "\t" + bp.getAffectedStatus());
		}
		fam.close();

		ArrayList<SNP> snpList = mapFile.getMarkerList();
		for (int i = 0; i < snpList.size(); i++)
		{
			SNP snp = snpList.get(i);
			bim.println(snp.getChromosome() + "\t" + snp.getName() + "\t" + snp.getDistance() + "\t" + snp.getPosition() + "\t" + snp.getFirstAllele() + "\t" + snp.getSecAllele());
		}
		bim.close();

		try
		{
			bedout.writeByte(byte1);
			bedout.writeByte(byte2);
			bedout.writeByte(byte3);

			for (int i = 0; i < gm.getNumMarker(); i++)
			{
				byte gbyte = 0;
				int idx = 0;

				for (int j = 0; j < gm.getNumIndivdial(); j++)
				{

					byte g = gm.getOriginalGenotypeScore(j, i);
					g <<= 2 * idx;
					gbyte |= g;
					idx++;

					if (j != (gm.getNumIndivdial() - 1))
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
			e.printStackTrace();
		}
	}
}
