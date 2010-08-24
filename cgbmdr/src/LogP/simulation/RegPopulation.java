package LogP.simulation;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;

import org.apache.commons.math.distribution.ChiSquaredDistribution;
import org.apache.commons.math.distribution.ChiSquaredDistributionImpl;

import LogP.LogPConstant;

import im.population.simulation.*;
import publicAccess.PublicData;
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
			populationSize = 200;
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
			rep = 100;
			if (param.size() > 15) {
				rep = Integer.parseInt(param.get(15));
			}
			// permutation
			permutation = 10;
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

	public static void main(String[] args) throws IOException {
		
		Parameter1 Param1 = null;
		if (args.length > 0) {
			Param1 = new Parameter1(args[0]);
		} else {
			String file = null;
			Param1 = new Parameter1(file);
		}

		/////////////////////////////////////////////////////////////Map
		MapMaker mapMaker = new MapMaker();
		double map[][] = { {0.0, 0.05392271280288696, 0.1356232762336731, 0.19633138179779053, 0.5293703675270081, 0.5840020775794983, 0.6356768608093262, 0.7160268425941467, 0.7968748807907104, 0.8703390955924988, 1.0 } };
		double[] len = {1.0};
		int[] m = {11};
		double[] criterion = {0.05, 0.25};
//		mapMaker.MakeRandomMap(len, m, Param1.seed, criterion);
//		map = mapMaker.getMap();
		
		/////////////////////////////////////////////////////////////QTL
		ArrayList QTL = new ArrayList();

		int[] chr0 = { 0 };
		double[] location0 = { 0.05 };
		int[] genotype0 = { 2 };
		double[] effect0 = { 0.5 };
		int environment1 = 0;
		AbstractLoci al0 = new AbstractLoci(chr0, location0, genotype0, effect0,
				environment1);
		QTL.add(al0);

		int[] chr1 = { 0 };
		double[] location1 = {0.45};
		int[] genotype1 = { 2 };
		double[] effect1 = { 0.5 };
		AbstractLoci al1 = new AbstractLoci(chr1, location1, genotype1, effect1,
				environment1);
		QTL.add(al1);

		
		
		double[] env = { 0.0 };
		double[][] weight = {{1, 0}};

		PrintStream Pout1 = null;
		ChiSquaredDistribution chi = new ChiSquaredDistributionImpl(weight.length);
		try {
			Pout1 = new PrintStream(new BufferedOutputStream(
					new FileOutputStream("LogAP.txt")));
		} catch (Exception E) {
			E.printStackTrace(System.err);
		}		
		
		ArrayList ControlMarker = null;
//		ControlMarker = mapMaker.getRandomMarker(map, Param1.seed);
		boolean isSelectMarker = true;
		int marker_control_stratagy = LogPConstant.SelectPairMarker;
		
		for (int i = 0; i < 100; i++) {
			AbstractMapping IM1 = new CompositeIntervalMapping(marker_control_stratagy);
			mapMaker.MakeRandomMap(len, m, Param1.seed + i, criterion);
			map = mapMaker.getMap();
			for (int j = 0; j < map.length; j++) {
				for (int k = 0; k < map[j].length; k++) {
					Pout1.print(map[i][j] + "\t");
				}
			}

			Param1.seed = i * 500;
			Param1.rep = 1;
			IM1.Simulation(Param1, QTL, env, weight, map, ControlMarker, isSelectMarker);
			IM1.SummuarySimulation();
		}
		Pout1.close();
	}
}
