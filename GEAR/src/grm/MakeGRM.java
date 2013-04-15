package grm;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;

import family.pedigree.Hukou;
import family.pedigree.file.MapFile;
import family.pedigree.file.PedigreeFile;
import family.pedigree.genotype.BPerson;
import family.plink.PLINKBinaryParser;
import family.plink.PLINKParser;
import family.popstat.GenotypeMatrix;
import family.qc.rowqc.SampleFilter;
import gear.Parameter;
import gear.util.FileProcessor;
import gear.util.Logger;
import sumstat.qc.rowqc.SumStatQC;

public class MakeGRM {
	private double maf_threshold = 1e-8;

	private GenotypeMatrix G;
	private int numMarker;
	private double[][] allelefreq;
	private double[][] genotypefreq;
	private MapFile snpMap;
	private SumStatQC ssQC;
	private PedigreeFile pf;
	
	public MakeGRM() {
		PLINKParser pp = null;
		if (Parameter.INSTANCE.getFileParameter().isSet()) {
			pp = new PLINKParser (Parameter.INSTANCE.getFileParameter().getPedFile(),
					              Parameter.INSTANCE.getFileParameter().getMapFile());
		} else if (Parameter.INSTANCE.getBfileParameter(0).isSet()) {
			pp = new PLINKBinaryParser (Parameter.INSTANCE.getBfileParameter(0).getBedFile(),
					                    Parameter.INSTANCE.getBfileParameter(0).getBimFile(),
					                    Parameter.INSTANCE.getBfileParameter(0).getFamFile());
		} else {
			Logger.printUserError("No input files.");
			System.exit(1);
		}
		pp.Parse();
		pf = pp.getPedigreeData();
		SampleFilter sf = new SampleFilter(pp.getPedigreeData(), pp.getMapData());
		SumStatQC ssQC = new SumStatQC(pp.getPedigreeData(), pp.getMapData(), sf);
		snpMap = pp.getMapData();
		GenotypeMatrix gm = new GenotypeMatrix(ssQC.getSample());
		G = gm;
		numMarker = G.getNumMarker();
		allelefreq = new double[numMarker][3];
	}

	public void makeGeneticRelationshipScore() {
		prepareMAF();
		StringBuffer sb = new StringBuffer();
		sb.append(Parameter.INSTANCE.out);
		PrintStream grm = null;
		BufferedWriter grmGZ = null;
		if (Parameter.INSTANCE.makeGRMTXTFlag) {
			sb.append(".grm.txt");
			grm = FileProcessor.CreatePrintStream(sb.toString());
		} else if (Parameter.INSTANCE.makeGRMFlag) {
			sb.append(".grm.gz");
			grmGZ = FileProcessor.ZipFielWriter(sb.toString());
		}

		for (int i = 0; i < G.getGRow(); i++) {
			for (int j = 0; j <= i; j++) {
				double[] s = GRMScore(i,j);
				if (Parameter.INSTANCE.makeGRMTXTFlag) {
					grm.println((i+1) + "\t" + (j+1) + "\t" + s[0] + "\t" + s[1]);
				} else {
					try {
						grmGZ.append((i+1) + "\t" + (j+1) + "\t" + s[0] + "\t" + s[1]);
					} catch (IOException e) {
						Logger.handleException(e, "error in writing '" + sb.toString() + "' for " + (i+1) + " " + (j+1) + ".");
					}
				}
			}
		}
		if (Parameter.INSTANCE.makeGRMTXTFlag) {
			grm.close();
		} else {
			try {
				grmGZ.close();
			} catch (IOException e) {
				Logger.handleException(e, "error in closing '" + sb.toString() + "'.");
			}
		}
		Logger.printUserLog("Writing GRM scores into '" + sb.toString() + "'.");		
		StringBuffer sb_id = new StringBuffer();
		sb_id.append(Parameter.INSTANCE.out);
		PrintStream grm_id = null;
		sb_id.append(".grm.id");
		grm_id = FileProcessor.CreatePrintStream(sb_id.toString());

		ArrayList<Hukou> H = pf.getHukouBook();
		for (int i = 0; i < H.size(); i++) {
			Hukou h = H.get(i);
			grm_id.println(h.getFamilyID() + "\t" + h.getIndividualID());
		}
		grm_id.close();
		Logger.printUserLog("Writing individual information into '" + sb_id.toString() + "'.");
	}

	private double[] GRMScore(int idx1, int idx2) {
		double[] s = { 0, 0 };

		for (int i = 0; i < allelefreq.length; i++) {

			int g1 = G.getAdditiveScore(idx1, i);
			int g2 = G.getAdditiveScore(idx2, i);
			double m = allelefreq[i][1];
			if (g1 == BPerson.MissingGenotypeCode
					|| g2 == BPerson.MissingGenotypeCode) {
				continue;
			} else {
				if (m < maf_threshold) {//ignor too little maf
					continue;
				} else {
					s[0]++;
					if (m < 0.5) {
						s[1] += (g1 - 2*m) * (g2 - 2*m)/(2 * m * (1-m));
					} else {
						s[1] += ((2-g1) - 2*(1-m)) * ((2-g2) - 2*(1-m))/(2 * m * (1-m));						
					}
				}
			}
		}

		if (s[0] > 0) {
			s[1] /= s[0];
		}

		return s;
	}

	private void prepareMAF() {
		if (Parameter.INSTANCE.ref_freq != null) {
			getRefFreq();

		} else {
			calAlleleFrequency();
		}
	}

	private void getRefFreq() {
		Logger.printUserLog("Got reference allele frequency from the " + Parameter.INSTANCE.ref_freq + ".\n");		
	}

	private void calAlleleFrequency() {
		int[][] g = G.getG();
		for (int i = 0; i < g.length; i++) {
			for (int j = 0; j < numMarker; j++) {
				int[] c = G.getBiAlleleGenotype(i, j);
				allelefreq[j][c[0]]++;
				allelefreq[j][c[1]]++;
				int idx = G.getAdditiveScore(i, j);
			}
		}

		for (int i = 0; i < numMarker; i++) {
			double wa = allelefreq[i][0] + allelefreq[i][1];
			double a = allelefreq[i][0] + allelefreq[i][1] + allelefreq[i][2];
			if (wa > 0) {
				for (int j = 0; j < allelefreq[i].length - 1; j++) {
					allelefreq[i][j] /= wa;
				}
				allelefreq[i][2] /= a;
			} else {
				allelefreq[i][2] = 1;
			}
		}		
		Logger.printUserLog("Calculated the reference allele frequency.");
	}

}
