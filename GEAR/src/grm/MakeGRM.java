package grm;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;

import family.pedigree.Hukou;
import family.pedigree.file.MapFile;
import family.pedigree.file.PedigreeFile;
import family.pedigree.file.SNP;
import family.pedigree.genotype.BPerson;
import family.plink.PLINKBinaryParser;
import family.plink.PLINKParser;
import family.popstat.GenotypeMatrix;
import family.qc.rowqc.SampleFilter;
import gear.Parameter;
import gear.util.FileProcessor;
import gear.util.Logger;
import gear.util.NewIt;
import gear.util.structure.MAF;
import sumstat.qc.rowqc.SumStatQC;

public class MakeGRM {
	private double maf_threshold = 1e-7;

	private GenotypeMatrix G;
	private int numMarker;
	private double[][] allelefreq;
	private MapFile snpMap;
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
						grmGZ.append((i+1) + "\t" + (j+1) + "\t" + s[0] + "\t" + s[1] + "\n");
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
				Logger.handleException(e, " error in closing '" + sb.toString() + "'.");
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

	public void makeGeneticRelationshipScore(int n0, int n1) {
		if (n0 <1 || n1 > (G.getGRow() * G.getGRow())/2) {
			Logger.printUserError("incorrect range for grm pairs : " + n0 + "," +n1);
			System.exit(0);
		} else {
			Logger.printUserLog("generating grm scores for the pairs from " + n0 + " to " + n1);
		}

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

		int i = 0, j = 0;
		for (; i < G.getGRow(); i++) {
			if( (i+1) * (i+1+1) / 2 >= n0) break;
		}
		j = n0 - i * (i+1)/2 - 1;
		int c = n0;

		while ( c <= n1 ) {
			double[] s = GRMScore(i,j);
			if (Parameter.INSTANCE.makeGRMTXTFlag) {
				grm.println((i+1) + "\t" + (j+1) + "\t" + s[0] + "\t" + s[1]);
			} else {
				try {
					grmGZ.append((i+1) + "\t" + (j+1) + "\t" + s[0] + "\t" + s[1] + "\n");
				} catch (IOException e) {
					Logger.handleException(e, "error in writing '" + sb.toString() + "' for " + (i+1) + " " + (j+1) + ".");
				}
			}
			if ( j < i) {
				j++;
			} else {
				i++;
				j = 0;
			}
			c++;
		}

		if (Parameter.INSTANCE.makeGRMTXTFlag) {
			grm.close();
		} else {
			try {
				grmGZ.close();
			} catch (IOException e) {
				Logger.handleException(e, " error in closing '" + sb.toString() + "'.");
			}
		}
		Logger.printUserLog("Writing GRM scores into '" + sb.toString() + "'.");		
		StringBuffer sb_id = new StringBuffer();
		sb_id.append(Parameter.INSTANCE.out);
		PrintStream grm_id = null;
		sb_id.append(".grm.id");
		grm_id = FileProcessor.CreatePrintStream(sb_id.toString());

		ArrayList<Hukou> H = pf.getHukouBook();
		for (int k = 0; k < H.size(); k++) {
			Hukou h = H.get(k);
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
		HashMap<String, MAF> refMap = NewIt.newHashMap();
		BufferedReader reader = FileProcessor.FileOpen(Parameter.INSTANCE.ref_freq);
		String line;
		try {
			line = reader.readLine();
			int idx = 0;
			while((line = reader.readLine())!=null) {
				line = line.trim();
				MAF maf = new MAF(line, ++idx);
				refMap.put(maf.getSNP(), maf);
			}
			Logger.printUserLog("Read " + idx + " SNPs in '" +  Parameter.INSTANCE.ref_freq + "'.\n");

		} catch (IOException e) {
			Logger.handleException(e, "An exception occurred when reading the maf file '" + Parameter.INSTANCE.ref_freq + "'.");
		}

		ArrayList<SNP> snpList = snpMap.getMarkerList();
		int c = 0;
		for (int i = 0; i < snpList.size(); i++) {
			SNP snp = snpList.get(i);
			String snpName = snp.getName();
			double f = 0;
			if(refMap.containsKey(snpName)) {
				f = refMap.get(snpName).getMAF();
			}
			if (allelefreq[i][1] <= 0.5) {
				if (allelefreq[i][1] < Parameter.INSTANCE.maf_range[0] || allelefreq[i][1] > Parameter.INSTANCE.maf_range[1]) {
					allelefreq[i][1] = 0;
				} else {
					allelefreq[i][1] = f;				
				}
			} else {
				if ( (1 - allelefreq[i][1]) < Parameter.INSTANCE.maf_range[0] || (1 - allelefreq[i][1]) > Parameter.INSTANCE.maf_range[1]) {
					allelefreq[i][1] = 0;
				} else {
					allelefreq[i][1] = f;				
				}
			}

			c++;
		}
		Logger.printUserLog("Got " + c + " matched reference alleles.\n");
	}

	private void calAlleleFrequency() {
		int[][] g = G.getG();
		for (int i = 0; i < g.length; i++) {
			for (int j = 0; j < numMarker; j++) {
				int[] c = G.getBiAlleleGenotype(i, j);
				allelefreq[j][c[0]]++;
				allelefreq[j][c[1]]++;
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
			if (allelefreq[i][1] <= 0.5) {
				if (allelefreq[i][1] < Parameter.INSTANCE.maf_range[0] || allelefreq[i][1] > Parameter.INSTANCE.maf_range[1]) {
					allelefreq[i][1] = 0;
				}
			} else {
				if ( (1 - allelefreq[i][1]) < Parameter.INSTANCE.maf_range[0] || (1 - allelefreq[i][1]) > Parameter.INSTANCE.maf_range[1]) {
					allelefreq[i][1] = 0;
				}				
			}
		}

		Logger.printUserLog("Calculated the reference allele frequency.");
	}

}
