package LogP.simulation;

import im.GenomeScan;
import im.IMBMatrix;
import im.IntervalPriorProbability;
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
import java.util.Set;

import publicAccess.ToolKit;
import regression.Likelihood;
import regression.LinearRegression;
import sun.security.provider.SystemSigner;
import LogP.PointMappingStatistic;
import LogP.Utils;
import LogP.simulation.RegPopulation.Parameter1;
import algorithm.CombinationGenerator;

import org.apache.commons.math.distribution.ChiSquaredDistribution;
import org.apache.commons.math.distribution.ChiSquaredDistributionImpl;

public class IntervalMapping extends AbstractMapping{

	public IntervalMapping() {
		super();
	}

	public void correctLogBon() {
		LogBon = -1 * Math.log10(0.05/Num_interval);
//		BonT = 0.05/1.574;
	}

	private void setup(ArrayList QTL, double[][]d, double[][] map) {
		int[][] pointIndex = new int[QTL.size()][3];
		pointPowerLOD = new Hashtable();
		pointPowerAT = new Hashtable();
		pointPowerDT = new Hashtable();
		pointPowerAPvalue = new Hashtable();
		pointPowerDPvalue = new Hashtable();
		for (int i = 0; i < QTL.size(); i++) {
			AbstractLoci al = (AbstractLoci) QTL.get(i);
			pointIndex[i][0] = al.getChr()[0];
			double location = d[al.getChr()[0]][al.getLoci()[0]];
			for (int j = 0; j < map[al.getChr()[0]].length-1; j++) {
				if(location < map[al.getChr()[0]][j+1]) {
					pointIndex[i][1] = j;
					pointIndex[i][2] = (new Double((location-map[al.getChr()[0]][j])*100)).intValue() + 1;
					break;
				}
			}
			pointPowerLOD.put(Utils.makeKey(pointIndex[i]), new Integer(0));
			pointPowerAT.put(Utils.makeKey(pointIndex[i]), new Integer(0));
			pointPowerDT.put(Utils.makeKey(pointIndex[i]), new Integer(0));
			pointPowerAPvalue.put(Utils.makeKey(pointIndex[i]), new Integer(0));
			pointPowerDPvalue.put(Utils.makeKey(pointIndex[i]), new Integer(0));
		}
	}

	protected void selectMarker() {

	}

	public ArrayList MappingProcedure(int i_rep, int isPermutation) {
		double[][] Y = new double[ap.IndividualNumber()][1];
		double[][] H0Matrix = null;
		H0Matrix = new double[ap.IndividualNumber()][1];
		for (int i = 0; i < H0Matrix.length; i++) {
			H0Matrix[i][0] = 1;
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
//		lmfull.MLE();
		ChiSquaredDistribution chi = new ChiSquaredDistributionImpl(weight.length);
		PrintStream Pout1 = null;
		try {
			Pout1 = new PrintStream(new BufferedOutputStream(
				new FileOutputStream("B.txt")));
		} catch (Exception E) {
			
		}
		for (int i = Param1.search_start; i <= Param1.search_end; i++) {
			imb.setOrder(i);
			CombinationGenerator cg = new CombinationGenerator(i, i, ap
					.SumIntevals());
			cg.generateCombination();
			List com = cg.get(i);

			for (Iterator e = com.iterator(); e.hasNext();) {
				String s = (String) e.next();
				int[] SNPIdx = ToolKit.StringToIntArray(s);
				double[][] Y_res = lmfull.getResponse();
				int[][] ChrInt = ChrInt(SNPIdx);
				double log0 = 0;

				LinearRegression H0lm = new LinearRegression(H0Matrix, Y_res);
				H0lm.MLE();
				Likelihood lkhd0 = new Likelihood(ap, gs, s);
				log0 = lkhd0.LogLikelihoodNullIM(H0lm);

				IntervalPriorProbability[] iip = new IntervalPriorProbability[SNPIdx.length];
				for (int j = 0; j < SNPIdx.length; j++) {
					iip[j] = gs.getIPPTable(ChrInt[j][0], ChrInt[j][1]);
				}
				for (int j = 0; j < iip.length; j++) {
					int steps = iip[j].getWalks();
					for (int jj = 0; jj < steps; jj++) {
						double[][] H1Matrix = imb.getICIMMatrixAtPoint(s,
								weight, jj);
						for(int k = 0; k < H1Matrix.length; k++) {
							Pout1.print(H1Matrix[k][1] + " ");
						}
						Pout1.println();
						LinearRegression H1lm = new LinearRegression(H1Matrix,
								Y_res);
						H1lm.MLE();
						Likelihood lkhd1 = new Likelihood(ap, gs, s);
						double log1 = lkhd1.LogLikelihoodAlternativeIM(H1lm,
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
							dominant_t = H1lm.getTStatic(2);
							dominant_t_p = H1lm.getPValueTTest(2);
						}
						double wald = H1lm.getWald(weight.length);
						double p_wald = 0;
						try {
							p_wald = chi.cumulativeProbability(wald);
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
		Pout1.close();
		return t;
	}
}
