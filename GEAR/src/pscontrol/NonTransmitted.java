package pscontrol;

import java.util.ArrayList;
import java.util.Iterator;

import parameter.AboutInfo;
import parameter.Parameter;
import pscontrol.hierarchy.AJHG2008;
import pscontrol.write.NonTransWriteBedSNPMajor;
import test.Test;
import util.NewIt;
import family.pedigree.PersonIndex;
import family.plink.PLINKBinaryParser;
import family.plink.PLINKParser;
import family.qc.rowqc.SampleFilter;

public class NonTransmitted {
	private String casualLociFile = null;
	private int[] casualLociIdx = null;
	private PLINKParser pp = null;
	private SampleFilter sf = null;
	private NonTransmitted QC = null;

	public NonTransmitted() {
		System.err.println("--nontrans procedure.");

//		if (Parameter.fileOption) {
//			pp = new PLINKParser(Parameter.pedfile, Parameter.mapfile);
//		}
		if (Parameter.INSTANCE.hasBFileOption()) {
			pp = new PLINKBinaryParser (Parameter.INSTANCE.getBedFile(),
					                    Parameter.INSTANCE.getBimFile(),
					                    Parameter.INSTANCE.getFamFile());
		} else {
			System.err.println("did not specify files.");
			Test.LOG.append("did not specify files.\n");
			Test.printLog();
			System.exit(0);
		}
		pp.Parse();
		sf = new SampleFilter(pp.getPedigreeData(), pp.getMapData());
	}

	public void GenerateNonTransmitted() {
		AJHG2008 ajhg2008 = new AJHG2008(pp.getPedigreeData(), pp.getMapData());
		ajhg2008.setSeed(Parameter.INSTANCE.nontransSeed);
		ajhg2008.RevvingUp(sf.getSample());

		ArrayList<PersonIndex> sample = ajhg2008.getSample();
		ArrayList<PersonIndex> ps = NewIt.newArrayList();

		System.err.println("Sample size: " + sample.size());
		for (Iterator<PersonIndex> e = sample.iterator(); e.hasNext();) {
			PersonIndex pi = e.next();
			if (pi.isPseudo()) {
				if (Parameter.INSTANCE.nontranscasesFlag) {
					if (pi.getPerson().getAffectedStatus().compareTo("2") != 0) {
						continue;
					}
				}
				if (Parameter.INSTANCE.nontranscontrolsFlag) {
					if (pi.getPerson().getAffectedStatus().compareTo("1") != 0) {
						continue;
					}
				}
				ps.add(pi);
			}
		}

		NonTransWriteBedSNPMajor writeSNP = new NonTransWriteBedSNPMajor(ps, ajhg2008
				.getMapFile().getMarkerList());
		StringBuilder out = new StringBuilder(); 
		out.append(Parameter.INSTANCE.out);
		out.append(".nt");
		writeSNP.WriteFile(out.toString());
	}
}
