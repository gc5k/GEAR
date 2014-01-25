package gear.pscontrol;

import java.util.ArrayList;
import java.util.Iterator;

import gear.CmdArgs;
import gear.family.pedigree.PersonIndex;
import gear.family.plink.PLINKBinaryParser;
import gear.family.plink.PLINKParser;
import gear.family.qc.rowqc.SampleFilter;
import gear.pscontrol.hierarchy.AJHG2008;
import gear.pscontrol.write.NonTransWriteBedSNPMajor;
import gear.util.Logger;
import gear.util.NewIt;

public class NonTransmitted
{

	private PLINKParser pp = null;
	private SampleFilter sf = null;

	public NonTransmitted()
	{
		Logger.printUserLog("--nontrans procedure.");

		if (CmdArgs.INSTANCE.getBFileArgs(0).isSet())
		{
			pp = new PLINKBinaryParser(CmdArgs.INSTANCE.getBFileArgs(0)
					.getBed(), CmdArgs.INSTANCE.getBFileArgs(0)
					.getBim(), CmdArgs.INSTANCE.getBFileArgs(0)
					.getFam());
		} else
		{
			Logger.printUserError("--bfile is not set.");
			System.exit(1);
		}
		pp.Parse();
		sf = new SampleFilter(pp.getPedigreeData(), pp.getMapData());
	}

	public void GenerateNonTransmitted()
	{
		AJHG2008 ajhg2008 = new AJHG2008(pp.getPedigreeData(), pp.getMapData());
		ajhg2008.setSeed(CmdArgs.INSTANCE.simuSeed);
		ajhg2008.RevvingUp(sf.getSample());

		ArrayList<PersonIndex> sample = ajhg2008.getSample();
		ArrayList<PersonIndex> ps = NewIt.newArrayList();

		Logger.printUserLog("Sample size: " + sample.size());
		for (Iterator<PersonIndex> e = sample.iterator(); e.hasNext();)
		{
			PersonIndex pi = e.next();
			if (pi.isPseudo())
			{
				if (CmdArgs.INSTANCE.nontranscasesFlag)
				{
					if (pi.getPerson().getAffectedStatus().compareTo("2") != 0)
					{
						continue;
					}
				}
				if (CmdArgs.INSTANCE.nontranscontrolsFlag)
				{
					if (pi.getPerson().getAffectedStatus().compareTo("1") != 0)
					{
						continue;
					}
				}
				ps.add(pi);
			}
		}

		NonTransWriteBedSNPMajor writeSNP = new NonTransWriteBedSNPMajor(ps,
				ajhg2008.getMapFile().getMarkerList());
		StringBuilder out = new StringBuilder();
		out.append(CmdArgs.INSTANCE.out);
		out.append(".nt");
		writeSNP.WriteFile(out.toString());
	}
}
