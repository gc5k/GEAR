package gear.pscontrol.write;

import java.io.DataOutputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Iterator;

import gear.CmdArgs;
import gear.ConstValues;
import gear.family.pedigree.PersonIndex;
import gear.family.pedigree.file.SNP;
import gear.family.pedigree.genotype.BPerson;
import gear.family.plink.PLINKBinaryParser;
import gear.family.plink.PLINKParser;
import gear.family.qc.rowqc.SampleFilter;
import gear.util.FileUtil;
import gear.util.Logger;

public class NonTransWriteBedSNPMajor
{

	private ArrayList<PersonIndex> PersonTable;
	private DataOutputStream os = null;
	private ArrayList<SNP> snpList;

	public NonTransWriteBedSNPMajor()
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

	public NonTransWriteBedSNPMajor(ArrayList<PersonIndex> pt, ArrayList<SNP> sl)
	{
		snpList = sl;
		PersonTable = pt;
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
		for (Iterator<PersonIndex> e = PersonTable.iterator(); e.hasNext();)
		{
			PersonIndex per = e.next();
			BPerson bp = per.getPerson();
			String[] id = bp.getPersonID().split("ajhg2008");
			pfam.append(bp.getFamilyID() + "\t" + id[0] + "\t" + bp.getDadID()
					+ "\t" + bp.getMomID() + "\t" + bp.getGender() + "\t"
					+ bp.getAffectedStatus() + "\n");
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
			Logger.handleException(e, "Cannot create file '" + sbed.toString()
					+ "'.");
		}

		try
		{
			os.writeByte(ConstValues.byte1);
			os.writeByte(ConstValues.byte2);
			os.writeByte(ConstValues.byte3);

			for (int i = 0; i < snpList.size(); i++)
			{
				byte gbyte = 0;
				int idx = 0;

				int posByte = i >> BPerson.shift;
				int posBit = (i & 0xf) << 1;

				for (int j = 0; j < PersonTable.size(); j++)
				{
					PersonIndex pi = PersonTable.get(j);
					BPerson bp = pi.getPerson();
					byte g = bp.getOriginalGenotypeScore(posByte, posBit);

					g <<= 2 * idx;
					gbyte |= g;
					idx++;

					if (j != (PersonTable.size() - 1))
					{
						if (idx == 4)
						{
							os.writeByte(gbyte);
							gbyte = 0;
							idx = 0;
						}
					} else
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
					"An exception occurred when writing file '"
							+ sbed.toString() + ".");
		}
	}

}
