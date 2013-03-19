package strand;

import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import parameter.Parameter;
import parameter.RegressionModel;
import parameter.AboutInfo;
import test.Test;
import util.FileProcessor;
import util.NewIt;
import family.pedigree.PersonIndex;
import family.pedigree.file.SNP;
import family.pedigree.genotype.BPerson;
import family.plink.PLINKBinaryParser;
import family.plink.PLINKParser;
import family.popstat.GenotypeMatrix;
import family.qc.rowqc.SampleFilter;

import util.SNPMatch;
import util.stat.Z;
import util.structure.Predictor;

public class MakePredictor {
	private GenotypeMatrix G1;

	private int[][] comSNPIdx;
	private double[][] allelefreq1;
	private double[] N1;
	private double[] N2;
	private ArrayList<Boolean> flag;

	private String[] title;
	private ArrayList<SNP> snpList1;
	private ArrayList<Predictor> predictorList = NewIt.newArrayList();

	ArrayList<Integer> scoreCoding = NewIt.newArrayList();

	private SampleFilter sf1;
	
	public MakePredictor() {
		System.err.print(AboutInfo.WELCOME_MESSAGE);
		
		readPredictor();

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
		snpList1 = sf1.getMapFile().getMarkerList();

	}

	public void BuildPredictor() {
		StringBuffer sb = new StringBuffer();
		sb.append(Parameter.INSTANCE.out);
		sb.append(".mergesnp");
		PrintStream ps = FileProcessor.CreatePrintStream(sb.toString());
		ps.append("SNP\tChr\tPos\tA1_1st\tA2_1st\tA1_2nd\tA2_2nd\tMAF_A1_1st\tMAF_A1_2nd\tFlip\tMerged\tP\tScheme\n");

		allelefreq1 = new double[G1.getNumMarker()][3];
		N1 = new double[G1.getNumMarker()];
		flag = NewIt.newArrayList();

		CalculateAlleleFrequency(G1, allelefreq1, N1);

		getCommonSNP(snpList1);

		int qualified_snp = 0;
		for (int i = 0; i < comSNPIdx[0].length; i++) {
			int scheme = 0;
			boolean ATGCLocus = false;
			boolean flip = false;

			SNP snp1 = snpList1.get(comSNPIdx[0][i]);
			Predictor maf2 = predictorList.get(comSNPIdx[1][i]);
			char a1_1 = snp1.getRefAllele();
			char a1_2 = snp1.getSecAllele();
			char a2_1 = maf2.getA1();

			boolean f = true;

			if (a1_1 == a2_1) {// scheme1
				scheme = 1;
				if (SNPMatch.Confusion(a1_1, a1_2)) {
					ATGCLocus = true;
				}
				scoreCoding.add(0);
			} else if (a1_2 == a2_1) {// scheme2
				scheme = 2;
				if (SNPMatch.Confusion(a1_1, a1_2)) {
					ATGCLocus = true;
				}
				scoreCoding.add(1);
			} else if (a1_1 == SNPMatch.Flip(a2_1)) {// scheme3
				scheme = 3;
				flip = true;
				scoreCoding.add(0);
			} else if (a1_2 == SNPMatch.Flip(a2_1)) {// scheme4
				scheme = 4;
				flip = true;
				scoreCoding.add(1);
			} else {// outlier
				scheme = 5;
				f = false;
				scoreCoding.add(0);
			}

			if (!Parameter.INSTANCE.keepATGC() && ATGCLocus) {
				f = false;
			}
			if (Parameter.INSTANCE.removeFlip() && flip) {
				f = false;
			}
			flag.add(f);
			if (f)
				qualified_snp++;

			ps.println(snpList1.get(comSNPIdx[0][i]).getName() + " "
					+ snpList1.get(comSNPIdx[0][i]).getChromosome() + " "
					+ snpList1.get(comSNPIdx[0][i]).getPosition() + " " + a1_1
					+ " " + a1_2 + " " + a2_1 + " " + " " + flip + " " + f
					+ " scheme" + scheme);
		}

		ps.close();
		if(qualified_snp == 0) {
			Test.LOG.append(qualified_snp + " common SNPs between two snp files.\nexit");
			System.err.println(qualified_snp + " common SNPs between two snp files.exit");
			System.exit(0);
		} else {
			Test.LOG.append(qualified_snp + " common SNPs can be used between two snp files.\n");
			System.err.println(qualified_snp + " common SNPs can be used between two snp files.");
		}

		WritePredictor();
	}

	private void getCommonSNP(ArrayList<SNP> snplist1) {
		HashMap<String, Integer> SNPMap = NewIt.newHashMap();
		for (Iterator<SNP> e = snplist1.iterator(); e.hasNext();) {
			SNP snp = e.next();
			SNPMap.put(snp.getName(), 0);
		}

		int c=0;
		HashMap<String, Integer> SNPMapList2 = NewIt.newHashMap();
		for (int i = 0; i < predictorList.size(); i++) {
			Predictor maf = predictorList.get(i);
			String snp_name = maf.getSNP();
			if (SNPMap.containsKey(snp_name)) {
				SNPMap.put(snp_name, 1);
				SNPMapList2.put(snp_name, i);
				c++;
			} else {
				SNPMap.put(snp_name, 0);
			}
		}

		if(c == 0) {
			Test.LOG.append(0 + " common SNPs between two snp files.\nexit");
			System.err.println(c + " common SNPs between two snp files.exit");
			System.exit(0);
		} else {
			Test.LOG.append(c + " common SNPs between two snp files.\n");
			System.err.println(c + " common SNPs between two snp files.");
		}

		comSNPIdx = new int[2][c];
		int idx1 = 0;
		for (int i = 0; i < snplist1.size(); i++ ) {
			SNP snp = snplist1.get(i);
			String snp_name = snp.getName();
			if (SNPMap.containsKey(snp_name) && SNPMap.get(snp_name).intValue() == 1) {
				comSNPIdx[0][idx1] = i;
				comSNPIdx[1][idx1] = SNPMapList2.get(snp_name).intValue();
				idx1++;
			}
		}
		System.out.println("idx1 "+ idx1);

	}

	public void CalculateAlleleFrequency(GenotypeMatrix G, double[][] frq, double[] n) {
		int[][] g = G.getG();
		for (int i = 0; i < g.length; i++) {
			for (int j = 0; j < G.getNumMarker(); j++) {
				int[] c = G.getBiAlleleGenotype(i, j);
				frq[j][c[0]]++;
				frq[j][c[1]]++;
			}
		}
		for (int i = 0; i < G.getNumMarker(); i++) {
			double w = frq[i][0] + frq[i][1];
			n[i] = frq[i][0] + frq[i][1] + frq[i][2];
			if (w > 0) {
				for (int j = 0; j < frq[i].length - 1; j++) {
					frq[i][j] /= w;
				}
				frq[i][2] /= n[i];
			} else {
				frq[i][2] = 1;
			}
		}
	}

	public void readPredictor() {

		BufferedReader reader = FileProcessor.FileOpen(Parameter.INSTANCE.getPredictorFile());
		String line;
		try {
			line = reader.readLine();
			int idx = 1;
			line = line.trim();
			title = line.split("\\s+");
			
			while((line = reader.readLine())!=null) {
				line = line.trim();
				Predictor maf = new Predictor(line, title.length, idx++);
				predictorList.add(maf);
			}
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	public void WritePredictor() {
		StringBuffer sbim = new StringBuffer();
		sbim.append(Parameter.INSTANCE.out);
		sbim.append(".predictor");
		PrintStream predictorFile = FileProcessor.CreatePrintStream(sbim.toString());
		
		int NMiss = 0;
		for (int i = 0; i < comSNPIdx[0].length; i++) {
			if(!flag.get(i)) {
				continue;
			}
			SNP snp = snpList1.get(comSNPIdx[0][i]);
			Predictor pd = predictorList.get(comSNPIdx[1][i]);
			if (Parameter.INSTANCE.isNA(pd.getField(Parameter.INSTANCE.getPredictorIdx()))) {
				NMiss++;
				continue;
			} else {
				double s = Double.parseDouble(pd.getField(Parameter.INSTANCE.getPredictorIdx()));
				if (Parameter.INSTANCE.getTranFunction() == RegressionModel.LINEAR) {
					if (scoreCoding.get(i).intValue() == 1) {
						s *= -1;
					}
				} else {
					if (s < 0) {
						NMiss++;
						continue;
					} else {
						if (scoreCoding.get(i).intValue() == 1) {
							s = Math.log(1/s);
						} else {
							s = Math.log(s);
						}
					}
				}
				predictorFile.append(snp.getName() + "\t" +  snp.getRefAllele() + "\t" + s +"\n");				
			}
		}
		predictorFile.close();
		System.out.println("write preditor to " + sbim.toString());
		System.out.println(NMiss + " snps have missing values and were not printed.");
	}
}
