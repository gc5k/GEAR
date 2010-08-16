package LogP.simulation;

import im.GenomeScan;
import im.IMBMatrix;
import im.IntervalPriorProbability;
import im.population.IMPopulation;
import im.population.simulation.AbstractPopulation;
import im.population.simulation.BackCrossPopulation;
import im.population.simulation.DHPopulation;
import im.population.simulation.F2Population;
import LogP.simulation.RegPopulation.Parameter1;

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import publicAccess.ToolKit;
import regression.Likelihood;
import regression.LinearRegression;

import algorithm.CombinationGenerator;

public class InclusiveIntervalMapping {
	
	ArrayList selectedMarker;
	int seed;
	AbstractPopulation ap;
	GenomeScan gs;
	Parameter1 Param1;
	int isPermutation;
	
	ArrayList threshold;
	ArrayList LOD;

	public InclusiveIntervalMapping(ArrayList sm, int s, GenomeScan g, Parameter1 p) {
		selectedMarker = sm;
		seed = s;
		gs = g;
		Param1 = p;
	}

	public ArrayList Simulation(ArrayList QTL, double[][] d, double[] env,
			int replication, boolean isLODTest) {
		for (int i_rep = 0; i_rep < replication; i_rep++) {
			if (Param1.permutation > 0 && i_rep == 0) {
				if (Param1.pt.compareTo("F2") == 0) {
					ap = new F2Population(Param1.populationSize,
							Param1.pheNum.length, Param1.pt, d, Param1.seed
									+ i_rep, Param1.mu, env, Param1.sd, QTL,
							Param1.mf);
				} else if (Param1.pt.compareTo("DH") == 0) {
					ap = new DHPopulation(Param1.populationSize,
							Param1.pheNum.length, Param1.pt, d, Param1.seed
									+ i_rep, Param1.mu, env, Param1.sd, QTL,
							Param1.mf);
				} else {
					ap = new BackCrossPopulation(Param1.populationSize,
							Param1.pheNum.length, Param1.pt, d, Param1.seed
									+ i_rep, Param1.mu, env, Param1.sd, QTL,
							Param1.mf);
				}
				ap.ProducePopulation();

				for (int i = 0; i < Param1.pheNum.length; i++) {
					ap.ProducePhenotype(i, Param1.MU, Param1.T);
				}

				threshold = new ArrayList();
				PrintStream Pout = null;
				PrintStream PoutCIM = null;
				try {
					Pout = new PrintStream(new BufferedOutputStream(
							new FileOutputStream("permuCIM.txt")));
					PoutCIM = new PrintStream(new BufferedOutputStream(
							new FileOutputStream("nullCIM.txt")));
				} catch (Exception E) {
					E.printStackTrace(System.err);
				}
				for (int i_permu = 0; i_permu < Param1.permutation; i_permu++) {
					ArrayList LOD = ICIM(i_rep, i_permu + 1);
					PoutCIM.print(LOD.size() + "\t");
					for (int ii = 0; ii < LOD.size(); ii++) {
						PoutCIM.print(LOD.get(ii) + ",");
					}
					PoutCIM.println();
					Double max = Collections.max(LOD);
					Double min = Collections.min(LOD);
					if (!isLODTest) {
						if (Math.abs(min.doubleValue()) > Math.abs(max
								.doubleValue())) {
							max = min;
						}
					}
					for (int ii = 0; ii < LOD.size(); ii++) {
						if (max == LOD.get(ii)) {
							Integer max_index = new Integer(ii);
						}
					}
					Pout.println(max);
					threshold.add(max);
				}
				Pout.close();
				PoutCIM.close();
				Collections.sort(threshold);
			}
			ArrayList L = ICIM(i_rep, 0);
			LOD.add(L);
		}
		return LOD;
	}

	public ArrayList ICIM(int i_rep, int isPermutation) {
		double[][] Y = new double[ap.IndividualNumber()][1];
		double[][] H0Matrix = null;
		H0Matrix = new double[ap.IndividualNumber()][1];
		for(int i = 0; i < H0Matrix.length; i++) {
			H0Matrix[i][0] = 1;
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
				double[][] Y_res = lmfull.quasiResidual(selectedMarker, SNPIdx[0]);
				int[][] ChrInt = ChrInt(ap, SNPIdx);
				double log0 = 0;
				LinearRegression H0lm = new LinearRegression(H0Matrix, Y_res);
				H0lm.MLE();
				Likelihood lkhd0 = new Likelihood(ap, gs, s);
				log0 = lkhd0.LogLikelihoodNullICIM(H0lm);
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

	public int[][] ChrInt(IMPopulation imp, int[] SNPIdx) {
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
}
