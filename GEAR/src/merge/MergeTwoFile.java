package merge;

import java.io.DataOutputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import parameter.Parameter;
import family.pedigree.PersonIndex;
import family.pedigree.file.SNP;
import family.pedigree.genotype.BPerson;
import family.plink.PLINKBinaryParser;
import family.plink.PLINKParser;
import family.popstat.GenotypeMatrix;
import family.qc.rowqc.SampleFilter;
import gear.util.FileProcessor;
import gear.util.Logger;
import gear.util.NewIt;
import gear.util.SNPMatch;
import gear.util.stat.Z;


public class MergeTwoFile {
	private GenotypeMatrix G1;
	private GenotypeMatrix G2;

	private int[][] comSNPIdx;
	private double[][] allelefreq1;
	private double[][] allelefreq2;
	private double[] N1;
	private double[] N2;
	private ArrayList<Boolean> flag;

	private ArrayList<SNP> snpList1;
	private ArrayList<SNP> snpList2;

	private ArrayList<PersonIndex> PersonTable1;
	private ArrayList<PersonIndex> PersonTable2;
	ArrayList<Integer> snpCoding = NewIt.newArrayList();

	private SampleFilter sf1;
	private SampleFilter sf2;

	
	private byte byte1 = 108;
	private byte byte2 = 27;
	private byte byte3 = 1;
	
	private DataOutputStream os = null;

	public MergeTwoFile() {
		PLINKParser pp1 = null;
		PLINKParser pp2 = null;
		if (Parameter.INSTANCE.hasBFileOption() && Parameter.INSTANCE.hasBFile2Option()) {
			pp1 = new PLINKBinaryParser (Parameter.INSTANCE.getBedFile(),
					                     Parameter.INSTANCE.getBimFile(),
					                     Parameter.INSTANCE.getFamFile());
			
			pp2 = new PLINKBinaryParser (Parameter.INSTANCE.getBedFile2(),
					                     Parameter.INSTANCE.getBimFile2(),
					                     Parameter.INSTANCE.getFamFile2());
		} else {
			Logger.printUserError("--bfile or --bfile2 is not set.");
			System.exit(1);
		}
		pp1.Parse();
		pp2.Parse();

		sf1 = new SampleFilter(pp1.getPedigreeData(), pp1.getMapData());
		sf2 = new SampleFilter(pp2.getPedigreeData(), pp2.getMapData());
		G1 = new GenotypeMatrix(sf1.getSample());
		G2 = new GenotypeMatrix(sf2.getSample());
		PersonTable1 = sf1.getSample();
		PersonTable2 = sf2.getSample();
		snpList1 = sf1.getMapFile().getMarkerList();
		snpList2 = sf2.getMapFile().getMarkerList();

	}

	public void Merge() {
		DecimalFormat fmt = new DecimalFormat("#.###E0");

		StringBuffer sb = new StringBuffer();
		sb.append(Parameter.INSTANCE.out);
		sb.append(".mergesnp");
		PrintStream ps = FileProcessor.CreatePrintStream(sb.toString());
		ps.append("SNP\tChr\tPos\tA1_1st\tA2_1st\tA1_2nd\tA2_2nd\tMAF_A1_1st\tMAF_A1_2nd\tFlip\tMerged\tP\tScheme\n");
		
		StringBuffer sb1 = new StringBuffer();
		sb1.append(Parameter.INSTANCE.out);
		sb1.append(".mergebadsnp");

		PrintStream ps1 = FileProcessor.CreatePrintStream(sb.toString());
		ps1.append("SNP\tChr\tPos\tA1_1st\tA2_1st\tA1_2nd\tA2_2nd\tMAF_A1_1st\tMAF_A1_2nd\tFlip\tMerged\tP\tScheme\n");

		allelefreq1 = new double[G1.getNumMarker()][3];
		N1 = new double[G1.getNumMarker()];
		allelefreq2 = new double[G2.getNumMarker()][3];
		N2 = new double[G2.getNumMarker()];
		flag = NewIt.newArrayList();

		CalculateAlleleFrequency(G1, allelefreq1, N1);
		CalculateAlleleFrequency(G2, allelefreq2, N2);

		getCommonSNP(snpList1, snpList2);

		int qualified_snp = 0;
		for (int i = 0; i < comSNPIdx[0].length; i++) {
			int scheme = 0;
			boolean ATGCLocus = false;
			boolean flip = false;
			SNP snp1 = snpList1.get(comSNPIdx[0][i]);
			SNP snp2 = snpList2.get(comSNPIdx[1][i]);
			char a1_1 = snp1.getRefAllele();
			char a1_2 = snp1.getSecAllele();
			char a2_1 = snp2.getRefAllele();
			char a2_2 = snp2.getSecAllele();

			double ref1 = allelefreq1[comSNPIdx[0][i]][0];
			double ref2 = allelefreq2[comSNPIdx[1][i]][0];
			boolean f = true;
			
			if (SNPMatch.IsBiallelic(a1_1, a1_2, a2_1, a2_2)) {
				if (a1_1 == a2_1) {//scheme1
					scheme = 1;
					if (SNPMatch.Confusion(a1_1, a1_2)) {
						ATGCLocus = true;
						if (ref1<0.5 && ref2<0.5) {
							if (ref1 < Parameter.INSTANCE.getMergeParameter().getMafCutoff() &&
								ref2 < Parameter.INSTANCE.getMergeParameter().getMafCutoff()) {
								f=true;
							} else {
								f=false;
							}
							flip = false;
							snpCoding.add(0);

						} else if (ref1<0.5 && ref2>0.5) {
							if (ref1 < Parameter.INSTANCE.getMergeParameter().getMafCutoff() &&
								ref2 > 1 - Parameter.INSTANCE.getMergeParameter().getMafCutoff()) {
								f=true;
							} else {
								f=false;
							}
							ref2 = 1-ref1;
							flip = true;
							snpCoding.add(1);

						} else if (ref1>0.5 && ref2<0.5) {
							if (ref1 > 1 - Parameter.INSTANCE.getMergeParameter().getMafCutoff() &&
								ref2 < Parameter.INSTANCE.getMergeParameter().getMafCutoff()) {
								f=true;
							} else {
								f=false;
							}
							ref2 = 1- ref2;
							flip = true;
							snpCoding.add(1);

						} else {
							if (ref1 > 1 - Parameter.INSTANCE.getMergeParameter().getMafCutoff() &&
								ref2 > 1 - Parameter.INSTANCE.getMergeParameter().getMafCutoff()) {
								f=true;
							} else {
								f=false;
							}
							flip = false;
							snpCoding.add(0);
						}
						//debug
						f = false;
					} else {
						flip = false;
						snpCoding.add(0);
						f=true;
					}

				} else if (a1_1 == a2_2) {//scheme2
					scheme = 2;
					if (SNPMatch.Confusion(a1_1, a1_2)) {
						ATGCLocus = true;
						if (ref1<0.5 && (1-ref2)<0.5) {
							if (ref1 < Parameter.INSTANCE.getMergeParameter().getMafCutoff() &&
								(1 - ref2) < Parameter.INSTANCE.getMergeParameter().getMafCutoff()) {
								f=true;
							} else {
								f=false;
							}
							snpCoding.add(1);
							ref2 = 1- ref2;
							flip = true;
						} else if (ref1<0.5 && (1-ref2) > 0.5) {
							if (ref1 < Parameter.INSTANCE.getMergeParameter().getMafCutoff() &&
								(1 - ref2) > 1 - Parameter.INSTANCE.getMergeParameter().getMafCutoff()) {
								f=true;
							} else {
								f=false;
							}
							flip = false;
							snpCoding.add(0);

						} else if (ref1>0.5 && (1-ref2) < 0.5) {
							if (ref1 > 1 - Parameter.INSTANCE.getMergeParameter().getMafCutoff() &&
								1 - ref2 < Parameter.INSTANCE.getMergeParameter().getMafCutoff()) {
								f=true;
							} else {
								f=false;
							}
							flip = false;
							snpCoding.add(0);

						} else {
							if (ref1 > 1 - Parameter.INSTANCE.getMergeParameter().getMafCutoff() &&
								1 - ref2 > 1 - Parameter.INSTANCE.getMergeParameter().getMafCutoff()) {
								f=true;
							} else {
								f=false;
							}
							flip = true;
							snpCoding.add(1);
							ref2 = 1- ref2;

						}
						//debug
						f = false;
					} else {
						flip = true;
						ref2 = 1- ref2;
						snpCoding.add(1);
						f=false;
					}
				} else if (a1_1 == SNPMatch.Flip(a2_1) ) {//scheme3
					scheme = 3;
					flip = true;
					snpCoding.add(0);
					f = true;
				} else if (a1_1 == SNPMatch.Flip(a2_2)) {//scheme4
					scheme = 4;
					flip = true;
					ref2 = 1-ref2;
					snpCoding.add(1);
					f=true;
					//debug
					f = false;
				} else {//outlier
					scheme = 5;
					f = false;
					snpCoding.add(0);
				}

				double p = Z.OddsRatioTestPvalueTwoTail(ref1, ref2, N1[comSNPIdx[0][i]], N2[comSNPIdx[1][i]]);
				if (p < Parameter.INSTANCE.getMergeParameter().getPCutoff()) {
					f = false;
				}
				if (!Parameter.INSTANCE.keepATGC() && ATGCLocus) {
					f = false;
				}
				if (Parameter.INSTANCE.removeFlip() && flip) {
					f = false;
				}
				flag.add(f);
				if (f) qualified_snp++;
				
				ps.println(snpList1.get(comSNPIdx[0][i]).getName() + " " + snpList1.get(comSNPIdx[0][i]).getChromosome() + " " + snpList1.get(comSNPIdx[0][i]).getPosition() + " " + a1_1 + " " + a1_2 + " "+ a2_1 + " " + a2_2 + " " + " " + fmt.format(ref1) + " " + fmt.format(allelefreq2[comSNPIdx[1][i]][0]) + " " + flip + " " + f + " " + p + " scheme" + scheme);
			} else {
				flag.add(false);
				snpCoding.add(0);
				ps1.println(snpList1.get(comSNPIdx[0][i]).getName() + " " + snpList1.get(comSNPIdx[0][i]).getChromosome() + " " + snpList1.get(comSNPIdx[0][i]).getPosition() + " " + a1_1 + " " + a1_2 + " " + a2_1 + " " + a2_2 + " " + " " + fmt.format(ref1) + " " + fmt.format(allelefreq2[comSNPIdx[1][i]][0]) + " " + flip + " " + f + " " + 1 + " scheme" + scheme);
			}
		}

		ps.close();
		ps1.close();
		if (qualified_snp == 0) {
			Logger.printUserError("Common SNPs between the two SNP files: None");
			System.exit(1);
		} else {
			Logger.printUserLog("Common SNP(s) between the two SNP files: " + qualified_snp);
		}

		WriteTwoFile();
	}

	private void getCommonSNP(ArrayList<SNP> snplist1, ArrayList<SNP> snplist2) {
		HashMap<String, Integer> SNPMap = NewIt.newHashMap();
		for (Iterator<SNP> e = snplist1.iterator(); e.hasNext();) {
			SNP snp = e.next();
			SNPMap.put(snp.getName(), 0);
		}

		int c=0;
		HashMap<String, Integer> SNPMapList2 = NewIt.newHashMap();
		for (int i = 0; i < snplist2.size(); i++) {
			SNP snp = snplist2.get(i);
			String snp_name = snp.getName();
			if (SNPMap.containsKey(snp_name)) {
				SNPMap.put(snp_name, 1);
				SNPMapList2.put(snp_name, i);
				c++;
			} else {
				SNPMap.put(snp_name, 0);
			}
		}

		if (c == 0) {
			Logger.printUserError("Common SNPs between the two SNP files: None");
			System.exit(1);
		} else {
			Logger.printUserLog("Common SNP(s) between the two SNP files:" + c);
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
			n[i] = frq[i][0] + frq[i][1];
			if (w > 0) {
				for (int j = 0; j < frq[i].length - 1; j++) {
					frq[i][j] /= w;
				}
				frq[i][2] /= frq[i][0] + frq[i][1] + frq[i][2];
			} else {
				frq[i][2] = 1;
			}
		}
	}

	
	public void WriteTwoFile() {
		StringBuffer sbim = new StringBuffer();
		sbim.append(Parameter.INSTANCE.out);
		sbim.append(".bim");
		PrintStream pbim = FileProcessor.CreatePrintStream(sbim.toString());
		
		for (int i = 0; i < comSNPIdx[0].length; i++) {
			if(!flag.get(i)) {
				continue;
			}
			SNP snp = snpList1.get(comSNPIdx[0][i]);
			pbim.append(snp.getChromosome() + "\t" +snp.getName() + "\t" + snp.getDistance() + "\t" + snp.getPosition() + "\t" + snp.getRefAllele() + "\t" + snp.getSecAllele() + "\n");
		}		
		pbim.close();
		
		StringBuffer sfam = new StringBuffer();
		sfam.append(Parameter.INSTANCE.out);
		sfam.append(".fam");
		PrintStream pfam = FileProcessor.CreatePrintStream(sfam.toString());		
		for (Iterator<PersonIndex> e = PersonTable1.iterator(); e.hasNext(); ) {
			PersonIndex per = e.next();
			BPerson bp = per.getPerson();
			pfam.append(bp.getFamilyID() + "\t" + bp.getPersonID() + "\t" + bp.getDadID() + "\t" + bp.getMomID() + "\t" + bp.getGender() + "\t" + bp.getAffectedStatus() + "\n");
		}

		for (Iterator<PersonIndex> e = PersonTable2.iterator(); e.hasNext(); ) {
			PersonIndex per = e.next();
			BPerson bp = per.getPerson();
			pfam.append(bp.getFamilyID() + "\t" + bp.getPersonID() + "\t" + bp.getDadID() + "\t" + bp.getMomID() + "\t" + bp.getGender() + "\t" + bp.getAffectedStatus() + "\n");
		}

		pfam.close();
		
		StringBuffer sbed = new StringBuffer();
		sbed.append(Parameter.INSTANCE.out);
		sbed.append(".bed");
		try {
			os = new DataOutputStream(new FileOutputStream(sbed.toString()));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		try {
			os.writeByte(byte1);
			os.writeByte(byte2);
			os.writeByte(byte3);

			for (int i = 0; i < comSNPIdx[0].length; i++) {
				if(!flag.get(i)) {
					continue;
				}

				int snpIdx = comSNPIdx[0][i];
				byte gbyte = 0;
				int idx = 0;

				int posByte = snpIdx >> BPerson.shift;
				int posBite = (snpIdx & 0xf) << 1;

				for (int j = 0; j < PersonTable1.size(); j++) {
					PersonIndex pi = PersonTable1.get(j);
					BPerson bp = pi.getPerson();
					byte g = bp.getOriginalGenotypeScore(posByte, posBite);

					g <<= 2 * idx;
					gbyte |= g;
					idx++;

					if (idx == 4) {
						os.writeByte(gbyte);
						gbyte = 0;
						idx = 0;
					}
				}

				int snpIdx2 = comSNPIdx[1][i];

				int posByte2 = snpIdx2 >> BPerson.shift;
				int posBite2 = (snpIdx2 & 0xf) << 1;

				for (int j = 0; j < PersonTable2.size(); j++) {
					PersonIndex pi = PersonTable2.get(j);
					BPerson bp = pi.getPerson();
					byte g = bp.getOriginalGenotypeScore(posByte2, posBite2);

					if(snpCoding.get(i).intValue() == 1) {
						switch(g) {
							case 0: g = 3; break;
							case 2: g = 2; break;
							case 3: g = 0; break;
							default: g = 1; break; //missing
						}
					}

					g <<= 2 * idx;
					gbyte |= (byte) g;
					idx++;

					if (j != (PersonTable2.size() - 1) ) {
						if (idx == 4) {
							os.writeByte(gbyte);
							gbyte = 0;
							idx = 0;
						}
					} else {
						os.writeByte(gbyte);
					}
				}
			}

			os.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
