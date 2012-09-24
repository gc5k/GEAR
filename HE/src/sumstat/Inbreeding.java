package sumstat;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;

import org.apache.commons.math.MathException;
import org.apache.commons.math.stat.inference.ChiSquareTestImpl;

import parameter.Parameter;
import sumstat.qc.rowqc.SumStatQC;
import test.Test;
import util.FileProcessor;
import util.NewIt;
import util.stat.FastFisherExactTest;

import family.pedigree.PersonIndex;
import family.pedigree.file.MapFile;
import family.pedigree.file.SNP;
import family.plink.PLINKBinaryParser;
import family.plink.PLINKParser;
import family.popstat.GenotypeMatrix;
import family.qc.rowqc.SampleFilter;

public class Inbreeding {
	private GenotypeMatrix G;
	private int numMarker;
	private double[] maf;
	private Parameter par;
	private MapFile snpMap;
	private SumStatQC ssQC;
	
	private double[] N;
	private String[][] indKeep;
	private HashMap<String, Integer> groupID = NewIt.newHashMap();
	private double[][] mafGroup;
	private int[] IndGroup;
	private double[][] w;
	private double[] Fst;

	private ArrayList<String> GroupInfor = NewIt.newArrayList();
	
	public Inbreeding(Parameter p) {
		par = p;

		PLINKParser pp = null;
		if (Parameter.fileOption) {
			pp = new PLINKParser(Parameter.pedfile, Parameter.mapfile);
		} else if (Parameter.bfileOption) {
			pp = new PLINKBinaryParser(Parameter.bedfile, Parameter.bimfile, Parameter.famfile);
		} else {
			System.err.println("did not specify files.");
			Test.LOG.append("did not specify files.\n");
			Test.printLog();
			System.exit(0);
		}
		pp.Parse();
		SampleFilter sf = new SampleFilter(pp.getPedigreeData(), pp.getMapData());
		ssQC = new SumStatQC(pp.getPedigreeData(), pp.getMapData(), sf);
		snpMap = pp.getMapData();
		GenotypeMatrix gm = new GenotypeMatrix(ssQC.getSample());
		setup(gm);
	}

	public void setup(GenotypeMatrix gm) {
		G = gm;
		numMarker = G.getNumMarker();
		readKeepFile();
		readGroupInfor();
	}

	public void CalculateFst() {
		int[][] g = G.getG();
		maf = new double[numMarker];
		mafGroup = new double[numMarker][groupID.size()];
		w = new double[numMarker][groupID.size()];
		N = new double[numMarker];
		
		Fst = new double[numMarker];
		for (int i = 0; i < g.length; i++) {
			int idx = IndGroup[i];
			for (int j = 0; j < numMarker; j++) {
				int[] c = G.getBiAlleleGenotype(i, j);
				if (c[0]==0) {
					maf[j]++;
					mafGroup[j][idx]++;
				}
				if (c[1]==0) {
					maf[j]++;
					mafGroup[j][idx]++;
				}
				if (c[1]!=2) {
					w[j][idx]++;
					N[j]++;
				}
			}
		}

		
		PrintWriter fstOut = null;
		try {
			fstOut = new PrintWriter(new String(par.out + ".fst"));
		} catch (IOException e) {
			e.printStackTrace();
		}
		fstOut.print("Chr\tSNP\tPos\tRefA"+ "\t");
		for (int i = 0; i < groupID.size(); i++) {
			fstOut.print("prop" + (i+1) + "\t" + "Freq" + (i+1) + "\t" + "NInd" + (i+1) + "\t");
		}
		fstOut.println("Freq\tFst");

		for (int i = 0; i < numMarker; i++) {
			double f = 0;
			maf[i] = maf[i]/(2*N[i]);
			for (int j = 0; j < groupID.size(); j++) {
				mafGroup[i][j] = mafGroup[i][j]/(2*w[i][j]);
				w[i][j] = w[i][j]/(N[i]);
			}

			fstOut.print(snpMap.getSNP(i).getChromosome()+ "\t" + snpMap.getMarkerName(i) + "\t" + snpMap.getSNP(i).getPosition()+ "\t" + snpMap.getSNP(i).getRefAllele()+ "\t");
			for (int j = 0; j < groupID.size(); j++) {
				f += w[i][j]*2*mafGroup[i][j]*(1-mafGroup[i][j]);
				fstOut.print(w[i][j] + "\t" + mafGroup[i][j] + "\t" + (int) (w[i][j]*N[i]) + "\t");
			}
			fstOut.print(maf[i] + "\t");
			if(maf[i]!=0) {
				Fst[i] = 1 - f/(2*maf[i]*(1-maf[i]));
			}
			fstOut.println(Fst[i]);
		}
		fstOut.close();
	}

	public void readGroupInfor() {
		ArrayList<PersonIndex> pt = ssQC.getSample();
		IndGroup = new int[pt.size()];
		Arrays.fill(IndGroup, -1);
		int i = 0;
		for (Iterator<PersonIndex> e = pt.iterator(); e.hasNext();) {
			PersonIndex pi = e.next();
			for (int j = 0; j < indKeep[0].length; j++) {
				if (pi.getFamilyID().compareTo(indKeep[0][j]) == 0 && pi.getIndividualID().compareTo(indKeep[1][j])==0 ) {
					String g = indKeep[2][j];
					int idx = groupID.get(g).intValue();
					IndGroup[i]=idx;
				}
			}
			i++;
		}
	}

	private void readKeepFile() {
		BufferedReader reader = FileProcessor.FileOpen(Parameter.fst_file);
		String line = null;
		ArrayList<String> famList = NewIt.newArrayList();
		ArrayList<String> indList = NewIt.newArrayList();
		ArrayList<String> groupList = NewIt.newArrayList();
		int gID = 0;
		try {
			while ((line = reader.readLine()) != null) {
				line.trim();
				String[] l = line.split(Parameter.whitespace);
				if(l.length < 3) continue;
				famList.add(l[0]);
				indList.add(l[1]);
				groupList.add(l[2]);
				if (!groupID.containsKey(l[2])) {
					groupID.put(l[2], gID);
					gID++;
				}
			}
		} catch (IOException e) {
			e.printStackTrace(System.err);
			System.exit(0);
		}
		indKeep = new String[3][];
		indKeep[0] = (String[]) famList.toArray(new String[0]);
		indKeep[1] = (String[]) indList.toArray(new String[0]);
		indKeep[2] = (String[]) groupList.toArray(new String[0]);
	}

}
