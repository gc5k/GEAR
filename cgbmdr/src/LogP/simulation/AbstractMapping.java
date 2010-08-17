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
import java.util.Map;
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
	Hashtable pointPowerEmpiricalAPvalue;
	Hashtable pointPowerEmpiricalDPvalue;

	ArrayList thresholdLOD;
	ArrayList thresholdATPos;
	ArrayList thresholdATNeg;
	ArrayList thresholdDT;
	ArrayList thresholdPvalue;

	ArrayList SimulationResults;

	double[][] distance;
	double LogBon;
	double t_at_0975 = 2.7611;
	double t_at_0025 = -2.784;
	double BonT;
	double typeI_bt;
	double typeI_EmpT;
	double typeI_LogBon;
	double EmpiricalLogBon;
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

	private void setup(ArrayList QTL, double[][] map, Parameter1 p) {
		int[][] pointIndex = new int[QTL.size()][3];
		pointPowerLOD = new Hashtable();
		pointPowerAT = new Hashtable();
		pointPowerDT = new Hashtable();
		pointPowerAPvalue = new Hashtable();
		pointPowerDPvalue = new Hashtable();
		pointPowerEmpiricalAPvalue = new Hashtable();
		pointPowerEmpiricalDPvalue = new Hashtable();
		for (int i = 0; i < QTL.size(); i++) {
			AbstractLoci al = (AbstractLoci) QTL.get(i);
			pointIndex[i][0] = al.getChr()[0];
			double location = al.getLocation()[0];
			for (int j = 0; j < map[al.getChr()[0]].length-1; j++) {
				if(location < map[al.getChr()[0]][j+1]) {
					pointIndex[i][1] = j;
					pointIndex[i][2] = (new Double((location-map[al.getChr()[0]][j])/p.step)).intValue()+ 1;
					break;
				}
			}
			pointPowerLOD.put(Utils.makeKey(pointIndex[i]), new Integer(0));
			pointPowerAT.put(Utils.makeKey(pointIndex[i]), new Integer(0));
			pointPowerDT.put(Utils.makeKey(pointIndex[i]), new Integer(0));
			pointPowerAPvalue.put(Utils.makeKey(pointIndex[i]), new Integer(0));
			pointPowerDPvalue.put(Utils.makeKey(pointIndex[i]), new Integer(0));
			pointPowerEmpiricalAPvalue.put(Utils.makeKey(pointIndex[i]), new Integer(0));
			pointPowerEmpiricalDPvalue.put(Utils.makeKey(pointIndex[i]), new Integer(0));
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

	public abstract void correctLogBon();
	
	public void selectMarker(ArrayList QTL, double[][] map) {
		selectedMarker = new ArrayList();
		for (Iterator e = QTL.iterator(); e.hasNext(); ) {
			AbstractLoci al = (AbstractLoci) e.next();
			int chr = al.getChr()[0];
			double loc = al.getLocation()[0];
			int c = 0;
			for (int i = 0; i < map.length; i++) {
				for (int j = 1; j < map[i].length; j++) {
					if (i == chr && map[i][j] > loc && map[i][j-1] < loc) {
						ArrayList mp1 = new ArrayList();
						mp1.add(new Integer(i));
						mp1.add(new Integer(j-1));
						ArrayList mp2 = new ArrayList();
						mp2.add(new Integer(i));
						mp2.add(new Integer(j));
						selectedMarker.add(mp1);
						selectedMarker.add(mp2);
						break;
					}
				}
			}
		}
	}

	public void Simulation(Parameter1 p, ArrayList QTL,
			double[] env, double[][] w, double[][] map, boolean isSelectMarker) {
		if (isSelectMarker ) {
			selectMarker(QTL, map);
		}
		weight = new double[w.length][];
		for (int i = 0; i < w.length; i++) {
			weight[i] = new double[w[i].length];
			System.arraycopy(w[i], 0, weight[i], 0, weight[i].length);
		}
		for (int i = 0; i < map.length; i++) {
			Num_interval += map[i].length == 1 ? map[i].length : (map[i].length - 1);
		}
		correctLogBon();
		makeFullMap(QTL, map);
		setup(QTL, map, p);
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
				PrintStream Pout1 = null;
				try {
					Pout1 = new PrintStream(new BufferedOutputStream(
							new FileOutputStream("Ghost.txt")));
				} catch (Exception E) {
					E.printStackTrace(System.err);
				}				
				thresholdLOD = new ArrayList();
				thresholdATPos = new ArrayList();
				thresholdATNeg = new ArrayList();
				thresholdDT = new ArrayList();
				thresholdPvalue = new ArrayList();
				for (int i_permu = 0; i_permu < Param1.permutation; i_permu++) {
					ArrayList IMStatistic = MappingProcedure(i_rep, i_permu + 1);
					ArrayList IMLOD = new ArrayList();
					ArrayList IMAT = new ArrayList();
					ArrayList IMDT = new ArrayList();
					ArrayList IMPvalue = new ArrayList();

					int c = 0;
					int ct = 0;
					int bt = 0;
					for (int ii = 0; ii < IMStatistic.size(); ii++) {
						PointMappingStatistic pms = (PointMappingStatistic) IMStatistic.get(ii);
						IMLOD.add(new Double(pms.get_LOD()));
						IMAT.add(new Double(pms.get_tStatistic_additive()));
						IMDT.add(new Double(pms.get_tStatistic_dominant()));
						IMPvalue.add(new Double(pms.get_logP_additive()));
						if (pms.get_logP_additive() > LogBon ) {
							c++;
						}
						if (pms.get_P_additive() < BonT) {
							bt++;
						}
						if(weight.length >1) {
							IMPvalue.add(new Double(pms.get_logP_dominance()));
						}
					}
					if (c>0) {
						typeI_LogBon +=1;
						for (Iterator ep = IMPvalue.iterator(); ep.hasNext(); ) {
							Pout1.print(ep.next() + "\t");
						}
						Pout1.println();						
					}
					thresholdLOD.add(Collections.max(IMLOD));

					thresholdATPos.add(Collections.max(IMAT));
					thresholdATNeg.add(Collections.min(IMAT));

					thresholdDT.add(Collections.max(IMDT));
					thresholdDT.add(Collections.min(IMDT));
					
					thresholdPvalue.add(Collections.max(IMPvalue));
				}
				Pout1.close();
			}
			ArrayList IMLOD = MappingProcedure(i_rep, 0);
			SimulationResults.add(IMLOD);
		}
	}
	
	public abstract ArrayList MappingProcedure(int i_rep, int isPermutation);

	private void calculatePower() {
		Collections.sort(thresholdPvalue);
		EmpiricalLogBon = ((Double) thresholdPvalue.get((new Double ((thresholdPvalue.size()-1)*0.95)).intValue())).doubleValue();
		System.out.println("BonT: " + BonT);
		System.out.println("LogBon: " + LogBon);
		System.out.println("EmpLogBon: " + EmpiricalLogBon);
		Collections.sort(thresholdLOD);
		double threshold_LOD = ((Double) thresholdLOD.get((new Double(
				(thresholdLOD.size()-1) * 0.95)).intValue())).doubleValue();
		System.out.println("threshold_LOD: " + threshold_LOD);

		Collections.sort(thresholdATPos);
		double threshold_at_975 = ((Double) thresholdATPos.get((new Double(
				(thresholdATPos.size()-1) * 0.975)).intValue())).doubleValue();
		double threshold_at_025 = ((Double) thresholdATNeg.get((new Double(
				(thresholdATNeg.size()-1) * 0.025)).intValue())).doubleValue();
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
				if (pointPowerEmpiricalAPvalue.containsKey(key) && pms.get_logP_additive() > EmpiricalLogBon) {
					int p = ((Integer) pointPowerEmpiricalAPvalue.get(key)).intValue();
					p++;
					pointPowerEmpiricalAPvalue.put(key, new Integer(p));
				}
				if (weight.length > 1) {
					if (pointPowerDPvalue.containsKey(key) && pms.get_logP_dominance() > LogBon) {
						int p = ((Integer) pointPowerDPvalue.get(key)).intValue();
						p++;
						pointPowerDPvalue.put(key, new Integer(p));
					}
					if (pointPowerEmpiricalDPvalue.containsKey(key) && pms.get_logP_dominance() > EmpiricalLogBon) {
						int p = ((Integer) pointPowerDPvalue.get(key)).intValue();
						p++;
						pointPowerEmpiricalDPvalue.put(key, new Integer(p));
					}
				}
				if (pointPowerAT.containsKey(key)) {
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
		System.out.println("type I LogBon: " + typeI_LogBon/Param1.permutation);	
		System.out.println("QTL\tLOD\tBonAP\tEmpAP\tBonDP\tEmpDP\tAT\tDT");
		for (Iterator e = keys.iterator(); e.hasNext();) {
			String key = (String) e.next();
			System.out.print(key + "\t");
			System.out.print((Integer) pointPowerLOD.get(key) + "\t");
			System.out.print((Integer) pointPowerAPvalue.get(key) + "\t");
			System.out.print((Integer) pointPowerEmpiricalAPvalue.get(key) + "\t");
			if (weight.length > 1) {
				System.out.print(pointPowerDPvalue.get(key) + "\t");
				System.out.print(pointPowerEmpiricalDPvalue.get(key) + "\t");
			} else {
				System.out.print("--\t--\t");
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
