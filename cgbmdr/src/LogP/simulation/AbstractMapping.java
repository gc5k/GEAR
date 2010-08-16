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
import java.util.Set;

import publicAccess.ToolKit;
import regression.Likelihood;
import regression.LinearRegression;
import sun.security.provider.SystemSigner;
import LogP.PointMappingStatistic;
import LogP.Utils;
import LogP.simulation.RegPopulation.Parameter1;
import algorithm.CombinationGenerator;

public abstract class AbstractMapping {
	ArrayList selectedMarker;
	AbstractPopulation ap;
	GenomeScan gs;
	Parameter1 Param1;
	int isPermutation;
	double[][] weight;
	Hashtable pointPowerLOD;
	Hashtable pointPowerAT;
	Hashtable pointPowerDT;
	Hashtable pointPowerAPvalue;
	Hashtable pointPowerDPvalue;

	ArrayList thresholdLOD;
	ArrayList thresholdAT;
	ArrayList thresholdDT;

	ArrayList SimulationResults;

	double[][] distance;
	double LogBon;
	double Num_interval;
	public AbstractMapping() {
		selectedMarker = null;
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
	private void setup(ArrayList QTL, double[][] map) {
		int[][] pointIndex = new int[QTL.size()][3];
		pointPowerLOD = new Hashtable();
		pointPowerAT = new Hashtable();
		pointPowerDT = new Hashtable();
		pointPowerAPvalue = new Hashtable();
		pointPowerDPvalue = new Hashtable();
		for (int i = 0; i < QTL.size(); i++) {
			AbstractLoci al = (AbstractLoci) QTL.get(i);
			pointIndex[i][0] = al.getChr()[0];
			double location = al.getLocation()[0];
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
	
	private void makeFullMap(ArrayList QTL, double[][] map) {
		ArrayList chromosomes = new ArrayList();
		for (int i = 0; i < map.length; i++) {
			ArrayList chromosome = new ArrayList();
			for (int j = 0; j < map[i].length; j++) {
				chromosome.add(new Double(map[i][j]));
			}
			chromosomes.add(chromosome);
		}
		for(Iterator e = QTL.iterator(); e.hasNext(); ) {
			AbstractLoci al = (AbstractLoci) e.next();
			int chr = al.getChr()[0];
			double loc = al.getLocation()[0];
			ArrayList chromosome = (ArrayList) chromosomes.get(chr);
			for (int j = 0; j < chromosome.size(); j++) {
				double loc_map = ((Double) chromosome.get(j)).doubleValue();
				if (loc<loc_map) {
					chromosome.add(j, new Double(loc));
					al.setLoci(j);
					break;
				}
			}
		}
		distance = new double[chromosomes.size()][];
		for (int i = 0; i < chromosomes.size(); i++ ) {
			ArrayList chromosome = (ArrayList) chromosomes.get(i);
			distance[i] = new double[chromosome.size()];
			for (int j = 0; j < chromosome.size(); j++) {
				distance[i][j] = ((Double) chromosome.get(j)).doubleValue();
			}
		}
	}

	public void Simulation(Parameter1 p, ArrayList QTL,
			double[] env, double[][] w, double[][] map, ArrayList sm) {
		selectedMarker = sm;
		weight = new double[w.length][];
		for (int i = 0; i < w.length; i++) {
			weight[i] = new double[w[i].length];
			System.arraycopy(w[i], 0, weight[i], 0, weight[i].length);
		}
		for (int i = 0; i < map.length; i++) {
			Num_interval += map[i].length - 1;
		}
		LogBon = -1 * Math.log10(0.05 / Num_interval);
		makeFullMap(QTL, map);
		setup(QTL, map);
		Param1 = p;
		SimulationResults = new ArrayList();
		for (int i_rep = 0; i_rep < Param1.rep; i_rep++) {
			if (Param1.pt.compareTo("F2") == 0) {
				ap = new F2Population(Param1.populationSize,
						Param1.pheNum.length, Param1.pt, distance,
						Param1.seed + i_rep, Param1.mu, env, Param1.sd, QTL,
						Param1.mf);
			} else if (Param1.pt.compareTo("DH") == 0) {
				ap = new DHPopulation(Param1.populationSize,
						Param1.pheNum.length, Param1.pt, distance,
						Param1.seed + i_rep, Param1.mu, env, Param1.sd, QTL,
						Param1.mf);
			} else {
				ap = new BackCrossPopulation(Param1.populationSize,
						Param1.pheNum.length, Param1.pt, distance,
						Param1.seed + i_rep, Param1.mu, env, Param1.sd, QTL,
						Param1.mf);
			}
			ap.ProducePopulation();

			for (int i = 0; i < Param1.pheNum.length; i++) {
				ap.ProducePhenotype(i, Param1.MU, Param1.T);
			}
			gs = new GenomeScan(ap, Param1.step);
			gs.CalculateIPP();

			if (Param1.permutation > 0 && i_rep == 0) {
				thresholdLOD = new ArrayList();
				thresholdAT = new ArrayList();
				thresholdDT = new ArrayList();

				for (int i_permu = 0; i_permu < Param1.permutation; i_permu++) {
					ArrayList IMStatistic = MappingProcedure(i_rep, i_permu + 1);
					ArrayList IMLOD = new ArrayList();
					ArrayList IMAT = new ArrayList();
					ArrayList IMDT = new ArrayList();

					for (int ii = 0; ii < IMStatistic.size(); ii++) {
						PointMappingStatistic pms = (PointMappingStatistic) IMStatistic.get(ii);
						IMLOD.add(new Double(pms.get_LOD()));
						IMAT.add(new Double(pms.get_tStatistic_additive()));
						IMDT.add(new Double(pms.get_tStatistic_dominant()));
					}
					Double maxLOD = Collections.max(IMLOD);
					thresholdLOD.add(maxLOD);
					
					thresholdAT.add(Collections.max(IMAT));
					thresholdAT.add(Collections.min(IMAT));

					thresholdDT.add(Collections.max(IMDT));
					thresholdDT.add(Collections.min(IMDT));
				}
			}
			ArrayList IMLOD = MappingProcedure(i_rep, 0);
			SimulationResults.add(IMLOD);
		}
	}
	
	public abstract ArrayList MappingProcedure(int i_rep, int isPermutation);

	private void calculatePower() {
		System.out.println("logBon: " + LogBon);
		Collections.sort(thresholdLOD);
		double threshold_LOD = ((Double) thresholdLOD.get((new Double(
				thresholdLOD.size() * 0.95)).intValue())).doubleValue();
		System.out.println("threshold_LOD: " + threshold_LOD);

		Collections.sort(thresholdAT);
		double threshold_at_975 = ((Double) thresholdAT.get((new Double(
				(thresholdAT.size()-1) * 0.975)).intValue())).doubleValue();
		double threshold_at_025 = ((Double) thresholdAT.get((new Double(
				(thresholdAT.size()-1) * 0.025)).intValue())).doubleValue();
		System.out.println("threshold_at_975: " + threshold_at_975);
		System.out.println("threshold_at_025: " + threshold_at_025);

		Collections.sort(thresholdDT);
		double threshold_dt_975 = ((Double) thresholdDT.get((new Double(
				(thresholdDT.size()-1) * 0.975)).intValue())).doubleValue();
		double threshold_dt_025 = ((Double) thresholdDT.get((new Double(
				(thresholdDT.size()-1) * 0.025)).intValue())).doubleValue();
		System.out.println("threshold_dt_975: " + threshold_dt_975);
		System.out.println("threshold_dt_025: " + threshold_dt_025);

		for (Iterator e = SimulationResults.iterator(); e.hasNext();) {
			ArrayList result = (ArrayList) e.next();
			for (Iterator e1 = result.iterator(); e1.hasNext();) {
				PointMappingStatistic pms = (PointMappingStatistic) e1.next();
				String key = pms.getKey();
				if (Param1.permutation > 0 && pointPowerLOD.containsKey(key) && pms.get_LOD() > threshold_LOD) {
					int p = ((Integer) pointPowerLOD.get(key)).intValue();
					p++;
					pointPowerLOD.put(key, new Integer(p));
				}
				if (pointPowerAPvalue.containsKey(key) && pms.get_logP_additive() > LogBon) {
					int p = ((Integer) pointPowerAPvalue.get(key)).intValue();
					p++;
					pointPowerAPvalue.put(key, new Integer(p));
				}
				if (weight.length > 1) {
					if (pointPowerDPvalue.containsKey(key)) {
						Integer power = (Integer) pointPowerDPvalue.get(key);
						if (pms.get_logP_dominant() > LogBon) {
							int p = power.intValue();
							p++;
							pointPowerDPvalue.put(key, new Integer(p));
						}
					}
				}
				if (pointPowerAT.containsKey(key)) {
					Integer power = (Integer) pointPowerAT.get(key);
					int p = ((Integer) pointPowerAT.get(key)).intValue();
					if (pms.get_tStatistic_additive() > 0
							&& pms.get_tStatistic_additive() > threshold_at_975) {
						p++;
						pointPowerAT.put(key, new Integer(p));
					}
					if (pms.get_tStatistic_additive() <= 0
							&& pms.get_tStatistic_additive() < threshold_at_025) {
						p++;
						pointPowerAT.put(key, new Integer(p));
					}
				}
				if (weight.length > 1) {
					if (pointPowerDT.contains(key)) {
						Integer power = (Integer) pointPowerDT.get(key);
						int p = ((Integer) pointPowerDT.get(key)).intValue();
						if (pms.get_tStatistic_dominant() > 0
								&& pms.get_tStatistic_dominant() > threshold_dt_975) {
							p++;
							pointPowerDT.put(key, new Integer(p));
						}
						if (pms.get_tStatistic_dominant() <= 0
								&& pms.get_tStatistic_dominant() > threshold_dt_025) {
							p++;
							pointPowerDT.put(key, new Integer(p));
						}
					}
				}
			}
		}
		Set keys = pointPowerLOD.keySet();
		System.out.println("QTL\tLOD\tAP\tDP\tAT\tDT");
		for (Iterator e = keys.iterator(); e.hasNext();) {
			String key = (String) e.next();
			System.out.print(key + "\t");
			System.out.print((Integer) pointPowerLOD.get(key) + "\t");
			System.out.print((Integer) pointPowerAPvalue.get(key) + "\t");
			if (weight.length > 1) {
				System.out.print(pointPowerDPvalue.get(key) + "\t");
			} else {
				System.out.print("--\t");
			}
			System.out.print(pointPowerAT.get(key) + "\t");
			if (weight.length > 1) {
				System.out.print(pointPowerDT.get(key) + "\t");
			} else {
				System.out.print("--\n");
			}
		}
	}

	public void SummuarySimulation() {
		calculatePower();
		PrintStream Pout1 = null;
		try {
			Pout1 = new PrintStream(new BufferedOutputStream(
					new FileOutputStream("logAPIM.txt")));
		} catch (Exception E) {
			E.printStackTrace(System.err);
		}
		for (int i = 0; i < SimulationResults.size(); i++) {
			ArrayList PointStatistics = (ArrayList) SimulationResults.get(i);
			for (int j = 0; j < PointStatistics.size(); j++) {
				PointMappingStatistic pms = (PointMappingStatistic) PointStatistics
						.get(j);
				Pout1.print(pms.get_logP_additive() + " ");
			}
			Pout1.println();
		}
		Pout1.close();

		try {
			Pout1 = new PrintStream(new BufferedOutputStream(
					new FileOutputStream("LODIM.txt")));
		} catch (Exception E) {
			E.printStackTrace(System.err);
		}
		for (int i = 0; i < SimulationResults.size(); i++) {
			ArrayList PointStatistics = (ArrayList) SimulationResults.get(i);
			for (int j = 0; j < PointStatistics.size(); j++) {
				PointMappingStatistic pms = (PointMappingStatistic) PointStatistics
						.get(j);
				Pout1.print(pms.get_LOD() + " ");
			}
			Pout1.println();
		}
		Pout1.close();
	}	
}
