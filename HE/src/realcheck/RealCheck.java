package realcheck;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.math.random.RandomDataImpl;

import parameter.Parameter;
import test.Test;
import util.FileProcessor;
import family.pedigree.PersonIndex;
import family.pedigree.file.SNP;
import family.pedigree.genotype.BPerson;
import family.plink.PLINKBinaryParser;
import family.plink.PLINKParser;
import family.popstat.GenotypeMatrix;
import family.qc.rowqc.SampleFilter;

public class RealCheck {
	private GenotypeMatrix G1;
	private GenotypeMatrix G2;

	private Parameter par;
	private int[] markerIdx;
	private ArrayList<SNP> snpList;

	private ArrayList<PersonIndex> PersonTable1;
	private ArrayList<PersonIndex> PersonTable2;
	
	private SampleFilter sf1;
	private SampleFilter sf2;
	public RealCheck(Parameter p) {
		System.err.print(Parameter.version);
		par = p;

		PLINKParser pp1 = null;
		PLINKParser pp2 = null;
		if (Parameter.bfileOption && Parameter.bfileOption2) {
			pp1 = new PLINKBinaryParser(Parameter.bedfile, Parameter.bimfile, Parameter.famfile);
			pp2 = new PLINKBinaryParser(Parameter.bedfile2, Parameter.bimfile2, Parameter.famfile2);
		} else {
			System.err.println("did not specify files.");
			Test.LOG.append("did not specify files.\n");
			Test.printLog();
			System.exit(0);
		}
		pp1.Parse();
		pp2.Parse();
		sf1 = new SampleFilter(pp1.getPedigreeData(), pp1.getMapData());
		sf2 = new SampleFilter(pp2.getPedigreeData(), pp2.getMapData());
		G1 = new GenotypeMatrix(sf1.getSample());
		G2 = new GenotypeMatrix(sf2.getSample());
		PersonTable1 = sf1.getSample();
		PersonTable2 = sf2.getSample();
		snpList = pp1.getMapData().getMarkerList();
	}

	public void Check() {
		StringBuffer sb = new StringBuffer();
		sb.append(par.out);
		sb.append(".real");
		PrintStream ps = FileProcessor.CreatePrintStream(sb.toString());
		
		getRandomMarker();
		for (int i = 0; i < G1.getGRow(); i++) {
			for (int j = 0; j < G2.getGRow(); j++) {
				double[] s = similarityScore(i, j);
				if (par.realcheckThresholdFlag) {
					if(s[0] > par.realcheckThreshold) {
						PersonIndex ps1 = PersonTable1.get(i);
						PersonIndex ps2 = PersonTable2.get(j);
						ps.print(ps1.getFamilyID() + " " + ps1.getIndividualID() + " " + ps2.getFamilyID() + " " + ps2.getIndividualID() + " " + s[0] + " " + s[1] + "/" + markerIdx.length + "\n");
					}
				}
			}
		}
		ps.close();
		
	}

	private double[] similarityScore(int idx1, int idx2) {
		double[] s = {0,0};

		for (int i = 0; i < markerIdx.length; i++) {
			int g1 = G1.getAdditiveScore(idx1, markerIdx[i]);
			int g2 = G2.getAdditiveScore(idx2, markerIdx[i]);
			if (g1 == BPerson.MissingGenotypeCode || g2 == BPerson.MissingGenotypeCode) continue;
			if (g1 == g2) {
				s[0]++;
			}
			s[1]++;
		}

		if(s[1]>0) {
			s[0] = s[0]/s[1];
		}

		return s;
	}
	
	public void getRandomMarker() {
		markerIdx = new int[par.realcheckMarkerNumber];
		RandomDataImpl rd = new RandomDataImpl();
		rd.reSeed(par.seed);
		markerIdx = rd.nextPermutation(snpList.size(), par.realcheckMarkerNumber);
		Arrays.sort(markerIdx);
		StringBuffer sb = new StringBuffer();
		sb.append(par.out);
		sb.append(".realsnp");

		PrintStream ps = FileProcessor.CreatePrintStream(sb.toString());
		for (int i = 0; i < markerIdx.length; i++) {
			SNP snp = snpList.get(markerIdx[i]);
			ps.print(snp.getChromosome() + " " + snp.getName() + " " + snp.getDistance() + " " + snp.getPosition() + " " + snp.getRefAllele() + " " + snp.getSecAllele() + "\n");
		}
		ps.close();
	}

}

