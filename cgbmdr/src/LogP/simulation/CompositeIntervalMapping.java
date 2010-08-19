package LogP.simulation;

import im.GenomeScan;
import im.IMBMatrix;
import im.IntervalPriorProbability;
import im.population.IMPopulation;
import im.population.simulation.AbstractLoci;
import im.population.simulation.AbstractPopulation;
import im.population.simulation.BackCrossPopulation;
import im.population.simulation.DHPopulation;
import im.population.simulation.F2Population;

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.math.distribution.ChiSquaredDistribution;
import org.apache.commons.math.distribution.ChiSquaredDistributionImpl;

import publicAccess.ToolKit;
import regression.Likelihood;
import regression.LinearRegression;
import LogP.PointMappingStatistic;
import LogP.Utils;
import LogP.simulation.RegPopulation.Parameter1;
import algorithm.CombinationGenerator;

import org.apache.commons.math.distribution.ChiSquaredDistribution;
import org.apache.commons.math.distribution.ChiSquaredDistributionImpl;

public class CompositeIntervalMapping extends AbstractMapping{

	public CompositeIntervalMapping() {
		super();
	}

	public void correctLogBon() {
		LogBon = (-1) * Math.log10(0.05/(Num_interval));
		BonT = 0.05/Num_interval;		
	}

	public ArrayList MappingProcedure(int i_rep, int isPermutation) {
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
		ChiSquaredDistribution chi = new ChiSquaredDistributionImpl(weight.length);
		for (int i = Param1.search_start; i <= Param1.search_end; i++) {
			imb.setOrder(i);
			CombinationGenerator cg = new CombinationGenerator(i, i, ap
					.SumIntevals());
			cg.generateCombination();
			List com = cg.get(i);
			for (Iterator e = com.iterator(); e.hasNext();) {
				String s = (String) e.next();
				double[][] H0Matrix = imb.getNullCIMMatrix(selectedMarker, s);
				LinearRegression H0lm = new LinearRegression(H0Matrix, Y);
				H0lm.MLE();
				Likelihood lkhd0 = new Likelihood(ap, gs, s);
				double log0 = lkhd0.LogLikelihoodNullCIM(H0lm);

				int[] SNPIdx = ToolKit.StringToIntArray(s);
				int[][] ChrInt = ChrInt(SNPIdx);
				IntervalPriorProbability[] iip = new IntervalPriorProbability[SNPIdx.length];

				for (int j = 0; j < SNPIdx.length; j++) {
					iip[j] = gs.getIPPTable(ChrInt[j][0], ChrInt[j][1]);
				}
				for (int j = 0; j < iip.length; j++) {
					int steps = iip[j].getWalks();
					for (int jj = 0; jj < steps; jj++) {
						double[][] H1Matrix = imb.getCIMMatrixAtPoint(
								selectedMarker, s, weight, jj);
						LinearRegression H1lm = new LinearRegression(H1Matrix,
								Y);
						H1lm.MLE();
						Likelihood lkhd1 = new Likelihood(ap, gs, s);
						double log1 = lkhd1.LogLikelihoodAlternativeCIM(H1lm,
								jj, weight);
						double lod = (log0 - log1) * (-1);
						double additive = H1lm.getCoefficient(1);
						double additive_sd = H1lm.getSD(1);
						double additive_t = H1lm.getTStatic(1);
						double additive_t_p = H1lm.getPValueTTest(1);
						double df = H1lm.DFResidual();
						double dominant = 0;
						double dominant_sd = 0;
						double dominant_t = 0;
						double dominant_t_p = 0;
						if (weight.length > 1) {
							dominant = H1lm.getCoefficient(2);
							dominant_sd = H1lm.getSD(2);
							dominant_t = H1lm.getPValueTTest(2);
							dominant_t_p = H1lm.getTStatic(2);
						}
						double wald = H1lm.getWald(weight.length);
						double p_wald = 0;
						try {
							chi.cumulativeProbability(wald);
						} catch (Exception E) {
							E.printStackTrace(System.err);
						}
						PointMappingStatistic pms = new PointMappingStatistic(
								ChrInt[j][0], ChrInt[j][1], jj, lod, additive,
								additive_sd, additive_t, additive_t_p,
								dominant, dominant_sd, dominant_t,
								dominant_t_p, df, wald, p_wald, weight.length, H1lm.getMSE());
						t.add(pms);
					}
				}
			}
		}
		return t;
	}

	public int[][] ChrInt(int[] SNPIdx) {
		int[][] chrint = new int[SNPIdx.length][2];
		int c = 0;
		int idx = 0;
		for (int i = 0; i < ap.ChromosomeNumber(); i++) {
			for (int j = 0; j < ap.IntervalNumberAtChromosome(i); j++) {
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
