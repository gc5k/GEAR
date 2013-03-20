package realcheck;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;

import org.apache.commons.math.random.RandomDataImpl;

import parameter.AboutInfo;
import parameter.Parameter;
import test.Test;
import util.FileProcessor;
import util.NewIt;
import util.TextHelper;
import family.pedigree.PersonIndex;
import family.pedigree.file.SNP;
import family.pedigree.genotype.BPerson;
import family.plink.PLINKBinaryParser;
import family.plink.PLINKParser;
import family.popstat.GenotypeMatrix;
import family.qc.rowqc.SampleFilter;

public class RealCheckOne {
	private GenotypeMatrix G1;

	private Parameter par;
	private int[] markerIdx;

	private ArrayList<SNP> snpList;

	private ArrayList<PersonIndex> PersonTable1;

	private SampleFilter sf1;

	public RealCheckOne() {
		System.err.print(AboutInfo.WELCOME_MESSAGE);

		PLINKParser pp1 = null;
		if (Parameter.INSTANCE.hasBFileOption()) {
			pp1 = new PLINKBinaryParser (Parameter.INSTANCE.getBedFile(),
					                     Parameter.INSTANCE.getBimFile(),
					                     Parameter.INSTANCE.getFamFile());
		} else {
			System.err.println("did not specify files.");
			Test.LOG.append("did not specify files.\n");
			Test.printLog();
			System.exit(0);
		}
		pp1.Parse();

		sf1 = new SampleFilter(pp1.getPedigreeData(), pp1.getMapData());
		G1 = new GenotypeMatrix(sf1.getSample());
		PersonTable1 = sf1.getSample();
		snpList = pp1.getMapData().getMarkerList();

	}

	public void Check() {

		StringBuffer sb = new StringBuffer();
		sb.append(par.out);
		sb.append(".real");
		PrintStream ps = FileProcessor.CreatePrintStream(sb.toString());

		if (Parameter.INSTANCE.getRealCheckParameter().getSnps() != null) {
			Test.LOG.append("generate similarity matrix with realcheck SNPs.\n");
			System.err.println("generate similarity matrix with realcheck SNPs.");
			getSelectedMarker();
		} else {
			getRandomMarker();
		}

		ps.print("file1.famid1 file1.id1 file1.famid2 file1.id2 score nmiss/loci\n");
		for (int i = 0; i < G1.getGRow(); i++) {
			for (int j = i; j < G1.getGRow(); j++) {
				double[] s = similarityScore(i, j);
				if (s[0] > Parameter.INSTANCE.getRealCheckParameter().getThresholdLower() && s[0] <= Parameter.INSTANCE.getRealCheckParameter().getThresholdUpper()) {
					PersonIndex ps1 = PersonTable1.get(i);
					PersonIndex ps2 = PersonTable1.get(j);
					ps.print(ps1.getFamilyID() + " "
							+ ps1.getIndividualID() + " "
							+ ps2.getFamilyID() + " "
							+ ps2.getIndividualID() + " " + s[0] + " "
							+ s[1] + "/" + markerIdx.length + "\n");
				}
			}
		}
		ps.close();

	}

	private double[] similarityScore(int idx1, int idx2) {
		double[] s = { 0, 0 };

		for (int i = 0; i < markerIdx.length; i++) {
			
			int idx = markerIdx[i];
			
			int g1 = G1.getAdditiveScore(idx1, idx);
			int g2 = G1.getAdditiveScore(idx2, idx);
			if (g1 == BPerson.MissingGenotypeCode
					|| g2 == BPerson.MissingGenotypeCode)
				continue;
			if (g1 == g2) {
				s[0]++;
			}
			s[1]++;
		}

		if (s[1] > 0) {
			s[0] = s[0] / s[1];
		}

		return s;
	}

	public void getSelectedMarker() {

		markerIdx = new int[sf1.getMapFile().getMarkerList().size()];

		for (int i = 0; i < markerIdx.length; i++) markerIdx[i] = i;

		StringBuffer sb = new StringBuffer();
		sb.append(par.out);
		sb.append(".realsnp");

		PrintStream ps = FileProcessor.CreatePrintStream(sb.toString());
		for (int i = 0; i < markerIdx.length; i++) {
			int idx = markerIdx[i];
			SNP snp = snpList.get(idx);
			ps.print(snp.getChromosome() + " " + snp.getName() + " "
					+ snp.getDistance() + " " + snp.getPosition() + " "
					+ snp.getRefAllele() + " " + snp.getSecAllele() + "\n");
		}
		ps.close();
	}

	public void getRandomMarker() {
		int mn = 0;
		int nMarker = sf1.getMapFile().getMarkerList().size();
		if (Parameter.INSTANCE.getRealCheckParameter().getMarkerNumber() > nMarker) {
			Test.LOG.append("realcheck marker number was reduced to " + nMarker + "\n");
			System.err.println("realcheck marker number was reduced to " + nMarker + "\n");
			mn = nMarker;
		} else {
			mn = Parameter.INSTANCE.getRealCheckParameter().getMarkerNumber();
		}

		markerIdx = new int[mn];
		RandomDataImpl rd = new RandomDataImpl();
		rd.reSeed(par.seed);

		markerIdx = rd.nextPermutation(markerIdx.length,
				mn);

		Arrays.sort(markerIdx);
		StringBuffer sb = new StringBuffer();
		sb.append(par.out);
		sb.append(".realsnp");

		PrintStream ps = FileProcessor.CreatePrintStream(sb.toString());
		for (int i = 0; i < markerIdx.length; i++) {
			int idx = markerIdx[i];
			SNP snp = snpList.get(idx);
			ps.print(snp.getChromosome() + " " + snp.getName() + " "
					+ snp.getDistance() + " " + snp.getPosition() + " "
					+ snp.getRefAllele() + " " + snp.getSecAllele() + "\n");
		}
		ps.close();
	}

	private ArrayList<String> readRealcheckSNPs() {
		BufferedReader reader = FileProcessor.FileOpen(Parameter.INSTANCE.getRealCheckParameter().getSnps());
		String line = null;
		ArrayList<String> selectedSNP = NewIt.newArrayList();
		try {
			while ((line = reader.readLine()) != null) {
				String[] l = line.split(TextHelper.WHITESPACE_DELIMITER);
				for (int i = 0; i < l.length; i++) {
					selectedSNP.add(l[i]);
				}
			}
		} catch (IOException e) {
			e.printStackTrace(System.err);
			System.exit(0);
		}
		Test.LOG.append("read " + selectedSNP.size() + " markers from " + Parameter.INSTANCE.getRealCheckParameter().getSnps() + ".\n");
		System.err.println("read " + selectedSNP.size() + " markers from " + Parameter.INSTANCE.getRealCheckParameter().getSnps() + ".");

		return selectedSNP;
	}
}
