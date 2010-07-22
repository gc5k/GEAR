import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.PrintStream;
import java.io.Reader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import org.apache.commons.math.linear.RealMatrix;

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
			switch2permutation = false;
			if (param.size() > 14) {
				switch2permutation = Boolean.parseBoolean(param.get(14));
			}
			// replication
			rep = 1;
			if (param.size() > 15) {
				rep = Integer.parseInt(param.get(15));
			}
			// permutation
			permutation = 100;
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
			int[] genotype1 = { 1 };
			double[] effect1 = { 0.5 };
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
		double d[][] = { { 0, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.65, 0.67,
				0.7, 0.75, 0.77, 0.8, 0.9, 1.0 } };
		Param2.ReadMap(d);

		// QTL
		ArrayList QTL = new ArrayList();

		int[] chr1 = { 0 };
		int[] loci1 = { 3 };
		int[] genotype1 = { 2 };
		double[] effect1 = { 1 };
		int environment1 = 0;
		AbstractLoci al1 = new AbstractLoci(chr1, loci1, genotype1, effect1,
				environment1);
		QTL.add(al1);

		int[] chr2 = { 0 };
		int[] loci2 = { 9 };
		int[] genotype2 = { 2 };
		double[] effect2 = { 1 };
		AbstractLoci al2 = new AbstractLoci(chr2, loci2, genotype2, effect2,
				environment1);
		QTL.add(al2);

		int[] chr3 = { 0 };
		int[] loci3 = { 11 };
		int[] genotype3 = { 2 };
		double[] effect3 = { -1 };
		AbstractLoci al3 = new AbstractLoci(chr3, loci3, genotype3, effect3,
				environment1);
		QTL.add(al3);

		Param2.ReadQTL(QTL);

		double[] env = { 0.0 };
		calculateMU(Param1, env, QTL, d);
		PrintStream Pout = new PrintStream(new BufferedOutputStream(
				new FileOutputStream("res1.txt")));
		for (int i_rep = 0; i_rep < Param1.rep; i_rep++) {
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

			GenomeScan gs = new GenomeScan(ap, Param1.step);
			gs.CalculateIPP();
			double[][] Y = new double[ap.IndividualNumber()][1];
			double[][] X;

			if (Param1.permutation == 0) {
				ArrayList ids = ap.getIDs();
				for (int i = 0; i < ids.size(); i++) {
					Integer id = (Integer) ids.get(i);
					Y[id.intValue()][0] = ap.PhenotypeAt(id.intValue(), 0);
				}
				IMBMatrix imb = new IMBMatrix(gs, ap);
				double[][] fm = imb.getFullMatrix();
				LinearRegression lm1 = new LinearRegression(fm, Y);
				lm1.MLE();
				for (int i = Param1.search_start; i <= Param1.search_end; i++) {
					imb.setOrder(i);
					CombinationGenerator cg = new CombinationGenerator(i, i, ap
							.SumIntevals());
					cg.generateCombination();
					List com = cg.get(i);
					double[][] Coeff = { { 1, 0 } };
					double[] mask = { 1, 0 };
					for (Iterator e = com.iterator(); e.hasNext();) {
						String s = (String) e.next();

						double[][] H0Matrix = imb.getNullCIMMatrix(s);
						LinearRegression H0lm = new LinearRegression(H0Matrix,
								Y);
						H0lm.MLE();
						Likelihood lkhd0 = new Likelihood(ap, gs, s);
						double log0 = lkhd0.LogLikelihood2(H0lm);

						int[] SNPIdx = ToolKit.StringToIntArray(s);
						int[][] ChrInt = ChrInt(ap, SNPIdx);
						IntervalPriorProbability[] iip = new IntervalPriorProbability[SNPIdx.length];
						for (int j = 0; j < SNPIdx.length; j++) {
							iip[j] = gs.getIPPTable(ChrInt[j][0], ChrInt[j][1]);
						}
						for (int j = 0; j < iip.length; j++) {
							int steps = iip[j].getWalks();
							for (int jj = 0; jj < steps; jj++) {
								double[][] H1Matrix = imb.getCIMMatrixAtPoint(
										s, Coeff, jj);
								LinearRegression H1lm = new LinearRegression(
										H1Matrix, Y);
								H1lm.MLE();
								Likelihood lkhd1 = new Likelihood(ap, gs, s);
								double log1 = lkhd1.LogLikelihood1(H1lm, jj);
								Pout.print((log0 - log1) * (-1) + " ");
							}
						}
					}
				}
			} else {
				for (int i_permu = 0; i_permu < Param1.permutation; i_permu++) {
					ap.Swith2Permutation(Param1.switch2permutation, Param1.seed
							* (i_rep * 100) + i_permu);
					ArrayList ids = ap.getIDs();
					for (int i = 0; i < ids.size(); i++) {
						Integer id = (Integer) ids.get(i);
						Y[id.intValue()][0] = ap.PhenotypeAt(id.intValue(), 0);
					}					
					IMBMatrix imb = new IMBMatrix(gs, ap);
					double[][] fm = imb.getFullMatrix();
					LinearRegression lm1 = new LinearRegression(fm, Y);
					lm1.MLE();
					for (int i = Param1.search_start; i <= Param1.search_end; i++) {
						imb.setOrder(i);
						CombinationGenerator cg = new CombinationGenerator(i,
								i, ap.SumIntevals());
						cg.generateCombination();
						List com = cg.get(i);
						double[][] Coeff = { { 1, 0 } };
						double[] mask = { 1, 0 };
						for (Iterator e = com.iterator(); e.hasNext();) {
							String s = (String) e.next();

							double[][] H0Matrix = imb.getNullCIMMatrix(s);
							LinearRegression H0lm = new LinearRegression(
									H0Matrix, Y);
							H0lm.MLE();
							Likelihood lkhd0 = new Likelihood(ap, gs, s);
							double log0 = lkhd0.LogLikelihood2(H0lm);

							int[] SNPIdx = ToolKit.StringToIntArray(s);
							int[][] ChrInt = ChrInt(ap, SNPIdx);
							IntervalPriorProbability[] iip = new IntervalPriorProbability[SNPIdx.length];
							for (int j = 0; j < SNPIdx.length; j++) {
								iip[j] = gs.getIPPTable(ChrInt[j][0],
										ChrInt[j][1]);
							}
							for (int j = 0; j < iip.length; j++) {
								int steps = iip[j].getWalks();
								for (int jj = 0; jj < steps; jj++) {
									double[][] H1Matrix = imb
											.getCIMMatrixAtPoint(s, Coeff, jj);
									LinearRegression H1lm = new LinearRegression(H1Matrix, Y);
									H1lm.MLE();
									Likelihood lkhd1 = new Likelihood(ap, gs, s);
									double log1 = lkhd1.LogLikelihood1(H1lm, jj);
									Pout.print((log0 - log1) * (-1) + " ");
								}
							}
						}
					}
					Pout.println();	
				}
				Pout.close();
			}
		}
	}
}