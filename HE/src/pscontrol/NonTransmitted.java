package pscontrol;

import java.util.ArrayList;
import java.util.Iterator;

import parameter.Parameter;
import pscontrol.hierarchy.AJHG2008;
import test.Test;
import util.NewIt;
import write.WriteBedSNPMajor;
import family.pedigree.PersonIndex;
import family.plink.PLINKBinaryParser;
import family.plink.PLINKParser;
import family.qc.rowqc.SampleFilter;

public class NonTransmitted {
	private String casualLociFile = null;
	private int[] casualLociIdx = null;
	private PLINKParser pp = null;
	private SampleFilter sf = null;
	Parameter par;

	public NonTransmitted (Parameter p) {
		System.err.print(Parameter.version);
		par = p;

		if (Parameter.fileOption) {
			pp = new PLINKParser(Parameter.pedfile, Parameter.mapfile);
		}
		if (Parameter.bfileOption) {
			pp = new PLINKBinaryParser(Parameter.bedfile, Parameter.bimfile, Parameter.famfile);
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
		ajhg2008.RevvingUp();
		
		ArrayList<PersonIndex> sample = ajhg2008.getSample();
		ArrayList<PersonIndex> ps = NewIt.newArrayList();
		
		for (Iterator<PersonIndex> e = sample.iterator(); e.hasNext(); ) {
			PersonIndex pi = e.next();
			if(pi.isPseudo()) {
				ps.add(pi);
			}
		}
		
		WriteBedSNPMajor writeSNP = new WriteBedSNPMajor(ps, ajhg2008.getMapFile().getMarkerList());
		writeSNP.WriteFile();
	}
}
