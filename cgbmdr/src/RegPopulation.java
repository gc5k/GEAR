import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.PrintStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import algorithm.*;
import im.population.simulation.AbstractLoci;
import im.population.IMPopulation;
import im.IMBMatrix;
import im.IntervalPriorProbability;
import im.GenomeScan;
import im.population.simulation.*;
import publicAccess.PublicData;
import publicAccess.ToolKit;
import regression.LinearRegression;
import regression.Likelihood;
import org.apache.commons.math.linear.RealMatrix;

/**
 * 
 * @author Guo-Bo Chen
 */
public class RegPopulation {

	public static class Parameter1 {
		protected int populationSize;// population size
		protected int[] pheNum;// number of phenotypes
		protected String pt;// population type
		protected double mu;// mean
		protected double sd;// sd
		protected int mf;// map function
		protected double step;// step of inteval
		protected boolean ssci;//
		protected int interval;
		protected int seed;
		protected double MU;// selection threshold
		protected double T; // threshold
		protected int search_start;
		protected int search_end;
		protected boolean switch2permutation;
		protected int rep;
		protected int permutation;

		public Parameter1(String file) {
			File config = null;
			if (file != null) {
				config = new File(file);
			}
			BufferedReader Reader = null;
			if (config != null) {
				try {
					Reader = new BufferedReader(new FileReader(config));
				} catch (Exception E) {
					E.printStackTrace(System.err);
				}
			}
			ArrayList<String> param = new ArrayList();
			String Line;
			try {
				while ((Reader != null) && (Line = Reader.readLine()) != null) {
					if (Line.length() == 0 || Line.startsWith("#")) {
						continue;
					}
					param.add(Line);
				}
			} catch (Exception E) {
				E.printStackTrace(System.err);
			}

			// population size
			populationSize = 300;
			if (param.size() > 0) {
				populationSize = Integer.parseInt(param.get(0));
			}
			// phenotype number
			pheNum = new int[1];
			pheNum[0] = 0;
			double[] os = { 0 };
			if (param.size() > 1) {
				pheNum[0] = Integer.parseInt(param.get(1));
				os[0] = 0;
			}
			// population type
			pt = new String("B1");
			if (param.size() > 2) {
				pt = param.get(2);
			}
			// population mu
			mu = 10;
			if (param.size() > 3) {
				mu = Double.parseDouble(param.get(3));
			}
			// residual
			sd = 1;
			if (param.size() > 4) {
				sd = Double.parseDouble(param.get(4));
			}
			// mapping function
			mf = 1;
			if (param.size() > 5) {
				mf = Integer.parseInt(param.get(5));
			}
			// step
			step = 0.01;
			if (param.size() > 6) {
				step = Double.parseDouble(param.get(6));
			}
			// search the same chromosome?
			ssci = true;
			if (param.size() > 7) {
				ssci = Boolean.parseBoolean(param.get(7));
			}
			// interval
			interval = 5;
			if (param.size() > 8) {
				interval = Integer.parseInt(param.get(8));
			}
			// seed
			seed = 82;
			if (param.size() > 9) {
				seed += Long.parseLong(param.get(9));
			}
			// selective mu
			MU = 0;
			if (param.size() > 10) {
				MU = Double.parseDouble(param.get(10));
			}
			// selective Threshold
			T = 0;
			if (param.size() > 11) {
				T = Double.parseDouble(param.get(11));
			}
			// search start
			search_start = 1;
			if (param.size() > 12) {
				search_start = Integer.parseInt(param.get(12));
			}
			// search end
			search_end = 1;
			if (param.size() > 13) {
				search_end = Integer.parseInt(param.get(13));
			}
			// switch to permutation
			switch2permutation = true;
			if (param.size() > 14) {
				switch2permutation = Boolean.parseBoolean(param.get(14));
			}
			// replication
			rep = 200;
			if (param.size() > 15) {
				rep = Integer.parseInt(param.get(15));
			}
			// permutation
			permutation = 10000;
			if (param.size() > 16) {
				permutation = Integer.parseInt(param.get(16));
			}
		}
	}

	public static class Parameter2 {
		protected ArrayList<String> param;

		public Parameter2(String file) {
			param = new ArrayList<String>();
			BufferedReader Reader2 = null;
			if (file != null) {
				try {
					Reader2 = new BufferedReader(new FileReader(file));
				} catch (Exception E) {
					E.printStackTrace(System.err);
				}
			}
			String Line2;
			try {
				while ((Reader2 != null)
						&& (Line2 = Reader2.readLine()) != null) {
					if (Line2.length() == 0 || Line2.startsWith("#")) {
						continue;
					}
					param.add(Line2);
				}
			} catch (Exception E) {
				E.printStackTrace(System.err);
			}
		}

		public void ReadMap(double[][] d) {
			if (param.size() > 1) {
				d = new double[Integer.parseInt(param.get(0))][];
				for (int k = 0; k < Integer.parseInt(param.get(0)); k++) {
					String[] distance = param.get(1 + k).split("[,\\s]++");
					d[k] = new double[distance.length];
					for (int kk = 0; kk < distance.length; kk++) {
						d[k][kk] = Double.parseDouble(distance[kk]);
					}
				}
			}
		}

		public void ReadQTL(ArrayList qtl) {
			int[] chr1 = { 0 };
			int[] loci1 = { 3 };
			int[] genotype1 = { 2 };
			double[] effect1 = { 1 };
			int pl = 0;
			if (param.size() > 0) {
				pl = Integer.parseInt(param.get(0));
				return;
			}
			if (param.size() > (pl + 1)) {
				qtl.clear();
				int qtlnumber = Integer.parseInt(param.get(pl + 1));
				for (int k = 0; k < qtlnumber; k++) {
					String[] chr = param.get((pl + 2 + k * 5))
							.split("[,\\s]++");
					chr1 = new int[chr.length];
					for (int kk = 0; kk < chr.length; kk++) {
						chr1[kk] = Integer.parseInt(chr[kk]);
					}

					String[] loc = param.get((pl + 2 + k * 5) + 1).split(
							"[,\\s]++");
					loci1 = new int[loc.length];
					for (int kk = 0; kk < loc.length; kk++) {
						loci1[kk] = Integer.parseInt(loc[kk]);
					}

					String[] fun = param.get((pl + 2 + k * 5) + 2).split(
							"[,\\s]++");
					genotype1 = new int[fun.length];
					for (int kk = 0; kk < fun.length; kk++) {
						genotype1[kk] = Integer.parseInt(fun[kk]);
					}

					String[] eff = param.get((pl + 2 + k * 5) + 3).split(
							"[,\\s]++");
					effect1 = new double[eff.length];
					for (int kk = 0; kk < eff.length; kk++) {
						effect1[kk] = Double.parseDouble(eff[kk]);
					}
					int envi1 = new Integer((param.get(pl + 2 + k * 5) + 4));
					AbstractLoci al = new AbstractLoci(chr1, loci1, genotype1,
							effect1, envi1);
					qtl.add(al);
				}
			}
		}
	}

	public static void calculateMU(Parameter1 Param1, double[] env,
			ArrayList QTL, double[][] d) {
		if ((Param1.MU - 0) > PublicData.epsilon) {
			AbstractPopulation ap;
			if (Param1.pt.compareTo("F2") == 0) {
				ap = new F2Population(Param1.populationSize,
						Param1.pheNum.length, Param1.pt, d, Param1.seed
								+ Param1.rep, Param1.mu, env, Param1.sd, QTL,
						Param1.mf);
			} else if (Param1.pt.compareTo("DH") == 0) {
				ap = new DHPopulation(Param1.populationSize,
						Param1.pheNum.length, Param1.pt, d, Param1.seed
								+ Param1.rep, Param1.mu, env, Param1.sd, QTL,
						Param1.mf);
			} else {
				ap = new BackCrossPopulation(Param1.populationSize,
						Param1.pheNum.length, Param1.pt, d, Param1.seed
								+ Param1.rep, Param1.mu, env, Param1.sd, QTL,
						Param1.mf);
			}
			ap.ProducePopulation();

			for (int i = 0; i < Param1.pheNum.length; i++) {
				ap.ProducePhenotype(i, Param1.MU, Param1.T);
			}
			Param1.MU = ap.getMean(0);
		}
	}

	public static int[][] ChrInt(IMPopulation imp, int[] SNPIdx) {
		int[][] chrint = new int[SNPIdx.length][2];
		int c = 0;
		int idx = 0;
		for (int i = 0; i < imp.ChromosomeNumber(); i++) {
			for (int j = 0; j < imp.IntervalNumberAtChromosome(i); j++) {
				if (c == SNPIdx[idx]) {
					chrint[idx][0] = i;
					chrint[idx][1] = j;
					idx++;
					if (idx == SNPIdx.length) {
						break;
					}
				}
				c++;
			}
			if (idx == SNPIdx.length) {
				break;
			}
		}
		return chrint;
	}

	public static void main(String[] args) throws IOException {
		Parameter1 Param1 = null;
		if (args.length > 0) {
			Param1 = new Parameter1(args[0]);
		} else {
			String file = null;
			Param1 = new Parameter1(file);
		}

		Parameter2 Param2 = null;
		if (args.length > 1) {
			Param2 = new Parameter2(args[1]);
		} else {
			String file = null;
			Param2 = new Parameter2(file);
		}
		double d[][] = { { 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55,
			               0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0 } };
		Param2.ReadMap(d);

		// QTL
		ArrayList QTL = new ArrayList();

		int[] chr1 = { 0 };
		int[] loci1 = { 6 };
		int[] genotype1 = { 2 };
		double[] effect1 = { 0.0179 };
		int environment1 = 0;
		AbstractLoci al1 = new AbstractLoci(chr1, loci1, genotype1, effect1,
				environment1);
		QTL.add(al1);

//		int[] chr2 = { 0 };
//		int[] loci2 = { 5 };
//		int[] genotype2 = { 2 };
//		double[] effect2 = { 0.25 };
//		AbstractLoci al2 = new AbstractLoci(chr2, loci2, genotype2, effect2,
//				environment1);
//		QTL.add(al2);
//
//		int[] chr3 = { 0 };
//		int[] loci3 = { 10 };
//		int[] genotype3 = { 2 };
//		double[] effect3 = { 0.25 };
//		AbstractLoci al3 = new AbstractLoci(chr3, loci3, genotype3, effect3,
//				environment1);
//		QTL.add(al3);

		Param2.ReadQTL(QTL);

		int[] PointIndex = {(new Double(d[0][loci1[0]] * 100)).intValue()};
							//(new Double(d[0][loci2[0]] * 100)).intValue(), 
							//(new Double(d[0][loci3[0]] * 100)).intValue()}; 
		
		ArrayList selectedMarker = null;
		//new ArrayList();

//		selectedMarker.add(new Integer(2));
		//selectedMarker.add(new Integer(7));

		boolean isLODTest = false;
		double[] powerIM = new double[PointIndex.length];
		double[] powerCIM = new double[PointIndex.length];
		double[] powertCIM= new double[PointIndex.length];
		double[] powerPermutationtCIM = new double[PointIndex.length];
		HashMap IM_permutation_Max_Index = new HashMap();
		HashMap CIM_permutation_Max_Index = new HashMap();
		HashMap tCIM_permutation_Max_Index = new HashMap();
		ArrayList<Double> thresholdIM = new ArrayList();
		ArrayList<Double> thresholdCIM = new ArrayList();
		ArrayList<Double> thresholdtCIM = new ArrayList();

		ArrayList<ArrayList<Double>> LODIM = new ArrayList();
		ArrayList<ArrayList<Double>> LODCIM = new ArrayList();
		ArrayList<ArrayList<Double>> LODtCIM = new ArrayList();

		double[] env = { 0.0 };
		calculateMU(Param1, env, QTL, d);

		for (int i_rep = 0; i_rep < Param1.rep; i_rep++) {
			//generate a population
			AbstractPopulation ap;
			if (Param1.pt.compareTo("F2") == 0) {
				ap = new F2Population(Param1.populationSize,
						Param1.pheNum.length, Param1.pt, d,
						Param1.seed + i_rep, Param1.mu, env, Param1.sd, QTL,
						Param1.mf);
			} else if (Param1.pt.compareTo("DH") == 0) {
				ap = new DHPopulation(Param1.populationSize,
						Param1.pheNum.length, Param1.pt, d,
						Param1.seed + i_rep, Param1.mu, env, Param1.sd, QTL,
						Param1.mf);
			} else {
				ap = new BackCrossPopulation(Param1.populationSize,
						Param1.pheNum.length, Param1.pt, d,
						Param1.seed + i_rep, Param1.mu, env, Param1.sd, QTL,
						Param1.mf);
			}
			ap.ProducePopulation();

			for (int i = 0; i < Param1.pheNum.length; i++) {
				ap.ProducePhenotype(i, Param1.MU, Param1.T);
			}
			PrintStream CIMfile = new PrintStream(new BufferedOutputStream(
					new FileOutputStream("CIMformat.mcd")));
			ap.CIMFormat(CIMfile);
			CIMfile.close();
			//config scan
			GenomeScan gs = new GenomeScan(ap, Param1.step);
			gs.CalculateIPP();

			//IM begins
			if (Param1.permutation > 0 && i_rep == 0) {
				thresholdIM = new ArrayList();
				PrintStream Pout = new PrintStream(new BufferedOutputStream(
						new FileOutputStream("permuIM.txt")));
				PrintStream PoutIM = new PrintStream(new BufferedOutputStream(
						new FileOutputStream("nullIIM.txt")));
				for (int i_permu = 0; i_permu < Param1.permutation; i_permu++) {
					ArrayList IMLOD = IM(selectedMarker, i_rep, ap, gs, Param1, i_permu + 1, true);
					PoutIM.print(IMLOD.size() + "\t");
					for (int ii = 0; ii < IMLOD.size(); ii++) {
						PoutIM.print(IMLOD.get(ii) + ",");
					}
					PoutIM.println();
					Double max = Collections.max(IMLOD);
					for (int ii = 0; ii < IMLOD.size(); ii++) {
						if (max == IMLOD.get(ii)) {
							Integer max_index = new Integer(ii);
							if(IM_permutation_Max_Index.containsKey(max_index)) {
								int c = ((Integer) IM_permutation_Max_Index.get(max_index)).intValue();
								c++;
								IM_permutation_Max_Index.put(max_index, new Integer(++c));
							} else {
								IM_permutation_Max_Index.put(max_index, new Integer(1));
							}
						}
					}
					Pout.println(max);
					thresholdIM.add(max);
				}
				Pout.close();
				PoutIM.close();
				Collections.sort(thresholdIM);
			}
			ArrayList IMLOD = IM(selectedMarker, i_rep, ap, gs, Param1, 0, true);
			LODIM.add(IMLOD);
			//IM ends
			
			//CIM t
			isLODTest = false;
			if (Param1.permutation > 0 && i_rep == 0) {
				thresholdtCIM = new ArrayList();
				PrintStream Pout = new PrintStream(new BufferedOutputStream(
						new FileOutputStream("permutCIM.txt")));
				PrintStream PoutCIM = new PrintStream(new BufferedOutputStream(
						new FileOutputStream("nulltCIM.txt")));
				for (int i_permu = 0; i_permu < Param1.permutation; i_permu++) {
					ArrayList tLOD = CIM(selectedMarker, i_rep, ap, gs, Param1, i_permu + 1, isLODTest);
					PoutCIM.print(tLOD.size()+ "\t");
					for (int ii = 0; ii < tLOD.size(); ii++) {
						PoutCIM.print(tLOD.get(ii) + ",");
					}
					PoutCIM.println();
					Double max = Collections.max(tLOD);
					Double min = Collections.min(tLOD);
					if (!isLODTest) {
						if (Math.abs(min.doubleValue())>Math.abs(max.doubleValue())) {
							max = min;
						}
					}
					for (int ii = 0; ii < tLOD.size(); ii++) {
						if (max == tLOD.get(ii)) {
							Integer max_index = new Integer(ii);
							if(tCIM_permutation_Max_Index.containsKey(max_index)) {
								int c = ((Integer) tCIM_permutation_Max_Index.get(max_index)).intValue();
								c++;
								tCIM_permutation_Max_Index.put(max_index, new Integer(++c));
							} else {
								tCIM_permutation_Max_Index.put(max_index, new Integer(1));
							}
						}
					}
					Pout.println(max);
					thresholdtCIM.add(max);
				}
				Pout.close();
				PoutCIM.close();
				Collections.sort(thresholdtCIM);
			}
			ArrayList tLOD = CIM(selectedMarker, i_rep, ap, gs, Param1, 0, isLODTest);
			LODtCIM.add(tLOD);
			//CIM end t
			
			//CIM begins LOD
			isLODTest = true;
			if (Param1.permutation > 0 && i_rep == 0) {
				thresholdCIM = new ArrayList();
				PrintStream Pout = new PrintStream(new BufferedOutputStream(
						new FileOutputStream("permuCIM.txt")));
				PrintStream PoutCIM = new PrintStream(new BufferedOutputStream(
						new FileOutputStream("nullCIM.txt")));
				for (int i_permu = 0; i_permu < Param1.permutation; i_permu++) {
					ArrayList LOD = CIM(selectedMarker, i_rep, ap, gs, Param1, i_permu + 1, isLODTest);
					PoutCIM.print(LOD.size()+ "\t");
					for (int ii = 0; ii < LOD.size(); ii++) {
						PoutCIM.print(LOD.get(ii) + ",");
					}
					PoutCIM.println();
					Double max = Collections.max(LOD);
					Double min = Collections.min(LOD);
					if (!isLODTest) {
						if (Math.abs(min.doubleValue())>Math.abs(max.doubleValue())) {
							max = min;
						}
					}
					for (int ii = 0; ii < LOD.size(); ii++) {
						if (max == LOD.get(ii)) {
							Integer max_index = new Integer(ii);
							if(CIM_permutation_Max_Index.containsKey(max_index)) {
								int c = ((Integer) CIM_permutation_Max_Index.get(max_index)).intValue();
								c++;
								CIM_permutation_Max_Index.put(max_index, new Integer(++c));
							} else {
								CIM_permutation_Max_Index.put(max_index, new Integer(1));
							}
						}
					}
					Pout.println(max);
					thresholdCIM.add(max);
				}
				Pout.close();
				PoutCIM.close();
				Collections.sort(thresholdCIM);
			}
			ArrayList LOD = CIM(selectedMarker, i_rep, ap, gs, Param1, 0, isLODTest);
			LODCIM.add(LOD);
			//CIM ends
		}

		PrintStream Pout1 = new PrintStream(new BufferedOutputStream(
				new FileOutputStream("LODtCIM.txt")));
		for (int i = 0; i < LODtCIM.size(); i++) {
			ArrayList tLod = (ArrayList) LODtCIM.get(i);
			for (int j = 0; j < tLod.size(); j++) {
				Pout1.print(tLod.get(j) + " ");
			}
			Pout1.println();
			for (int j = 0; j < PointIndex.length; j++) {
				if(Param1.permutation == 0) {
					continue;
				}
				if ((Double) tLod.get(PointIndex[j])>0) {
					if(((Double) tLod.get(PointIndex[j])) > ((Double) thresholdtCIM.get((int)(0.975*Param1.permutation)))) {
						powerPermutationtCIM[j] += 1;
					}
					if(((Double) tLod.get(PointIndex[j])).doubleValue() > (1-0.05/(d[0].length))) {
						powertCIM[j] += 1;;
					}
				} else {
					if(((Double) tLod.get(PointIndex[j])) > ((Double) thresholdtCIM.get((int)(0.025*Param1.permutation)))) {
						powerPermutationtCIM[j] += 1;
					}
					if(((Double) tLod.get(PointIndex[j])).doubleValue() > (0.05/(d[0].length))) {
						powertCIM[j] += 1;;
					}
				}
			}
		}
		Pout1.close();
		PrintStream Pout = new PrintStream(new BufferedOutputStream(
				new FileOutputStream("LODCIM.txt")));
		for (int i = 0; i < LODCIM.size(); i++) {
			ArrayList Lod = (ArrayList) LODCIM.get(i);
			for (int j = 0; j < Lod.size(); j++) {
				Pout.print(Lod.get(j) + " ");
			}
			Pout.println();
			for (int j = 0; j < PointIndex.length; j++) {
				if(Param1.permutation == 0) {
					continue;
				}
				if(((Double) Lod.get(PointIndex[j])) > ((Double) thresholdCIM.get((int)(0.95*Param1.permutation)))) {
					powerCIM[j] += 1;
				}
			}
		}
		Pout.close();

		PrintStream Pout2 = new PrintStream(new BufferedOutputStream(
				new FileOutputStream("LODIM.txt")));
		for (int i = 0; i < LODIM.size(); i++) {
			ArrayList Lod = (ArrayList) LODIM.get(i);
			for (int j = 0; j < Lod.size(); j++) {
				Pout2.print(Lod.get(j) + " ");
			}
			Pout2.println();
			for (int j = 0; j < PointIndex.length; j++) {
				if(Param1.permutation == 0) {
					continue;
				}
				if(((Double) Lod.get(PointIndex[j])) > ((Double) thresholdIM.get((int)(0.95*Param1.permutation)))) {
					powerIM[j] += 1;
				}
			}
		}
		Pout2.close();
		
//		System.out.println( "Threshold ICIM: "+thresholdICIM.get((int)(0.95*Param1.permutation)));
		for (int i = 0; i < PointIndex.length; i++) {
			System.out.println(powertCIM[i] + " " + powerPermutationtCIM[i] + " " + powerCIM[i]);
		}
		System.out.println("++++++++++++++++");
		Set Ikeys = tCIM_permutation_Max_Index.keySet();
		for (Iterator e = Ikeys.iterator(); e.hasNext(); ) {
			Integer Idx = (Integer) e.next();
			System.out.println(Idx + "->" + tCIM_permutation_Max_Index.get(Idx));
		}
		System.out.println("++++++++++++++++");
		Set keys = CIM_permutation_Max_Index.keySet();
		for (Iterator e = keys.iterator(); e.hasNext(); ) {
			Integer Idx = (Integer) e.next();
			System.out.println(Idx + "->" + CIM_permutation_Max_Index.get(Idx));
		}
	}

	public static ArrayList ICIM(ArrayList selectedMarker, int i_rep, IMPopulation ap, GenomeScan gs, Parameter1 Param1, int isPermutation, boolean shouldKeepMean) {
		double[][] Y = new double[ap.IndividualNumber()][1];
		double[][] H0Matrix = null;
		if (shouldKeepMean) {
			H0Matrix = new double[ap.IndividualNumber()][1];
			for(int i = 0; i < H0Matrix.length; i++) {
				H0Matrix[i][0] = 1;
			}
		}
		ArrayList t = new ArrayList();
		ArrayList ids;
		if (isPermutation == 0) {
			ids = ap.getIDs();
		} else {
			ap.Swith2Permutation(Param1.switch2permutation, Param1.seed * 10000
					 + isPermutation );
			ids = ap.getPermutatedIDs();
		}
		for (int i = 0; i < ids.size(); i++) {
			Integer id = (Integer) ids.get(i);
			if (isPermutation == 0) {
				Y[id.intValue()][0] = ap.PhenotypeAt(id.intValue(), 0);
			} else {
				Y[i][0] = ap.PhenotypeAt(id.intValue(), 0);
			}
		}
		IMBMatrix imb = new IMBMatrix(gs, ap);
		double[][] fm = imb.getFullMatrix(selectedMarker);
		LinearRegression lmfull = new LinearRegression(fm, Y);
		lmfull.MLE();
		for (int i = Param1.search_start; i <= Param1.search_end; i++) {
			imb.setOrder(i);
			CombinationGenerator cg = new CombinationGenerator(i, i, ap.SumIntevals());
			cg.generateCombination();
			List com = cg.get(i);
			double[][] Coeff = { {1, 0} };
			for (Iterator e = com.iterator(); e.hasNext();) {
				String s = (String) e.next();
				int[] SNPIdx = ToolKit.StringToIntArray(s);
				double[][] Y_res = lmfull.quasiResidual(selectedMarker, SNPIdx[0], shouldKeepMean);
				int[][] ChrInt = ChrInt(ap, SNPIdx);
				double log0 = 0;
				if (shouldKeepMean) {
					LinearRegression H0lm = new LinearRegression(H0Matrix, Y_res);
					H0lm.MLE();
					Likelihood lkhd0 = new Likelihood(ap, gs, s);
					log0 = lkhd0.LogLikelihoodNullICIM(H0lm);
				}
				IntervalPriorProbability[] iip = new IntervalPriorProbability[SNPIdx.length];
				for (int j = 0; j < SNPIdx.length; j++) {
					iip[j] = gs.getIPPTable(ChrInt[j][0], ChrInt[j][1]);
				}
				for (int j = 0; j < iip.length; j++) {
					int steps = iip[j].getWalks();
					for (int jj = 0; jj < steps; jj++) {
						double[][] H1Matrix = imb.getICIMMatrixAtPoint(s, Coeff,
								jj);
						LinearRegression H1lm = new LinearRegression(H1Matrix,
								Y_res);
						H1lm.MLE();
						Likelihood lkhd1 = new Likelihood(ap, gs, s);
						double log1 = lkhd1.LogLikelihoodAlternativeICIM(H1lm, jj);
						t.add(new Double((log0 - log1) * (-1)));
					}
				}
			}
		}
		return t;
	}

	public static ArrayList CIM(ArrayList selectedMarker, int i_rep, IMPopulation ap, GenomeScan gs,
			Parameter1 Param1, int isPermutation, boolean isLODTest) {
		double[][] Y = new double[ap.IndividualNumber()][1];
		ArrayList ids;
		ArrayList t = new ArrayList();
		if (isPermutation == 0) {
			ids = ap.getIDs();
		} else {
			ap.Swith2Permutation(Param1.switch2permutation, Param1.seed
					* 100000 + isPermutation);
			ids = ap.getPermutatedIDs();
		}
		for (int i = 0; i < ids.size(); i++) {
			Integer id = (Integer) ids.get(i);
			if (isPermutation == 0) {
				Y[id.intValue()][0] = ap.PhenotypeAt(id.intValue(), 0);
			} else {
				Y[i][0] = ap.PhenotypeAt(id.intValue(), 0);
			}
		}
		IMBMatrix imb = new IMBMatrix(gs, ap);
		double[][] fm = imb.getFullMatrix(selectedMarker);
		LinearRegression lm1 = new LinearRegression(fm, Y);
		lm1.MLE();
		for (int i = Param1.search_start; i <= Param1.search_end; i++) {
			imb.setOrder(i);
			CombinationGenerator cg = new CombinationGenerator(i, i, ap
					.SumIntevals());
			cg.generateCombination();
			List com = cg.get(i);
			double[][] Coeff = { { 1, 0 } };
			for (Iterator e = com.iterator(); e.hasNext();) {
				String s = (String) e.next();
				double[][] H0Matrix = imb.getNullCIMMatrix(selectedMarker, s);
				LinearRegression H0lm = new LinearRegression(H0Matrix, Y);
				H0lm.MLE();
				Likelihood lkhd0 = new Likelihood(ap, gs, s);
				double log0 = lkhd0.LogLikelihoodNullCIM(H0lm);

				int[] SNPIdx = ToolKit.StringToIntArray(s);
				int[][] ChrInt = ChrInt(ap, SNPIdx);
				IntervalPriorProbability[] iip = new IntervalPriorProbability[SNPIdx.length];

				for (int j = 0; j < SNPIdx.length; j++) {
					iip[j] = gs.getIPPTable(ChrInt[j][0], ChrInt[j][1]);
				}
				for (int j = 0; j < iip.length; j++) {
					int steps = iip[j].getWalks();
					for (int jj = 0; jj < steps; jj++) {
						double[][] H1Matrix = imb.getCIMMatrixAtPoint(selectedMarker, s, Coeff,
								jj);
						LinearRegression H1lm = new LinearRegression(H1Matrix,
								Y);
						H1lm.MLE();
						Likelihood lkhd1 = new Likelihood(ap, gs, s);
						double log1 = lkhd1.LogLikelihoodAlternativeCIM(H1lm, jj);
						if (isLODTest) {
							t.add(new Double((log0 - log1) * (-1)));
						} else {
							t.add(H1lm.getPValueTTest(1));														
						}
					}
				}
			}
		}
		return t;
	}

	public static ArrayList IM(ArrayList selectedMarker, int i_rep, IMPopulation ap, GenomeScan gs, Parameter1 Param1, int isPermutation, boolean shouldKeepMean) {
		double[][] Y = new double[ap.IndividualNumber()][1];
		double[][] H0Matrix = null;
		if (shouldKeepMean) {
			H0Matrix = new double[ap.IndividualNumber()][1];
			for(int i = 0; i < H0Matrix.length; i++) {
				H0Matrix[i][0] = 1;
			}
		}
		ArrayList t = new ArrayList();
		ArrayList ids;
		if (isPermutation == 0) {
			ids = ap.getIDs();
		} else {
			ap.Swith2Permutation(Param1.switch2permutation, Param1.seed
					* (i_rep * 100) + isPermutation);
			ids = ap.getPermutatedIDs();
		}
		for (int i = 0; i < ids.size(); i++) {
			Integer id = (Integer) ids.get(i);
			if (isPermutation == 0) {
				Y[id.intValue()][0] = ap.PhenotypeAt(id.intValue(), 0);
			} else {
				Y[i][0] = ap.PhenotypeAt(id.intValue(), 0);
			}
		}
		IMBMatrix imb = new IMBMatrix(gs, ap);
		double[][] fm = imb.getFullMatrix(selectedMarker);
		LinearRegression lmfull = new LinearRegression(fm, Y);
		lmfull.MLE();
		for (int i = Param1.search_start; i <= Param1.search_end; i++) {
			imb.setOrder(i);
			CombinationGenerator cg = new CombinationGenerator(i, i, ap.SumIntevals());
			cg.generateCombination();
			List com = cg.get(i);
			double[][] Coeff = { {1, 0} };
			for (Iterator e = com.iterator(); e.hasNext();) {
				String s = (String) e.next();
				int[] SNPIdx = ToolKit.StringToIntArray(s);
				double[][] Y_res = lmfull.getResponse();
				int[][] ChrInt = ChrInt(ap, SNPIdx);
				double log0 = 0;
				if (shouldKeepMean) {
					LinearRegression H0lm = new LinearRegression(H0Matrix, Y_res);
					H0lm.MLE();
					Likelihood lkhd0 = new Likelihood(ap, gs, s);
					log0 = lkhd0.LogLikelihoodNullICIM(H0lm);
				}
				IntervalPriorProbability[] iip = new IntervalPriorProbability[SNPIdx.length];
				for (int j = 0; j < SNPIdx.length; j++) {
					iip[j] = gs.getIPPTable(ChrInt[j][0], ChrInt[j][1]);
				}
				for (int j = 0; j < iip.length; j++) {
					int steps = iip[j].getWalks();
					for (int jj = 0; jj < steps; jj++) {
						double[][] H1Matrix = imb.getICIMMatrixAtPoint(s, Coeff,
								jj);
						LinearRegression H1lm = new LinearRegression(H1Matrix,
								Y_res);
						H1lm.MLE();
						Likelihood lkhd1 = new Likelihood(ap, gs, s);
						double log1 = lkhd1.LogLikelihoodAlternativeICIM(H1lm, jj);
						t.add(new Double((log0 - log1) * (-1)));
					}
				}
			}
		}
		return t;
	}
	
}
