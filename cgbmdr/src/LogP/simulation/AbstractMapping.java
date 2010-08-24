package LogP.simulation;

import im.GenomeScan;
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
import java.util.Set;

import LogP.PointMappingStatistic;
import LogP.Utils;
import LogP.simulation.RegPopulation.Parameter1;

import org.apache.commons.math.distribution.ChiSquaredDistribution;
import org.apache.commons.math.distribution.ChiSquaredDistributionImpl;
import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math.stat.descriptive.DescriptiveStatisticsImpl;
import org.apache.commons.math.stat.regression.SimpleRegression;

public abstract class AbstractMapping {
	ArrayList selectedMarker;
	AbstractPopulation ap;
	GenomeScan gs;
	Parameter1 Param1;
	int isPermutation;
	double[][] weight;
	Hashtable pointPowerF;
	Hashtable pointPowerWald;
	Hashtable pointPowerLOD;
	Hashtable pointPowerAT;
	Hashtable pointPowerDT;
	Hashtable pointPowerAPvalue;
	Hashtable pointPowerDPvalue;
	Hashtable pointPowerEmpiricalAPvalue;
	Hashtable pointPowerEmpiricalDPvalue;

	ArrayList thresholdPF;
	ArrayList thresholdLOD;
	ArrayList thresholdAP;
	ArrayList thresholdDP;
	ArrayList thresholdPwald;
	ArrayList SimulationResults;

	double threshold_LOD;
	double threshold_lodp;
	double threshold_Pwald;
	double threshold_PF;
	
	double[][] distance;
	double LogBon;
	double t_at_0975 = 2.7611;
	double t_at_0025 = -2.784;
	double BonT;
	double typeI_bt;
	double typeI_EmpT;
	double typeI_LogBon;
	double EmpiricalLogBonAP;
	double EmpiricalLogBonDP;
	double Num_interval;

	ArrayList QTL;
	double[][] map;
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
		pointPowerF = new Hashtable();
		pointPowerWald = new Hashtable();
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
			for (int j = 0; j < map[al.getChr()[0]].length - 1; j++) {
				if (location < map[al.getChr()[0]][j + 1]) {
					pointIndex[i][1] = j;
					pointIndex[i][2] = (new Double(
							(location - map[al.getChr()[0]][j]) / p.step))
							.intValue() + 1;
					break;
				}
			}
			pointPowerF.put(Utils.makeKey(pointIndex[i]), new Integer(0));
			pointPowerWald.put(Utils.makeKey(pointIndex[i]), new Integer(0));
			pointPowerLOD.put(Utils.makeKey(pointIndex[i]), new Integer(0));
			pointPowerAT.put(Utils.makeKey(pointIndex[i]), new Integer(0));
			pointPowerDT.put(Utils.makeKey(pointIndex[i]), new Integer(0));
			pointPowerAPvalue.put(Utils.makeKey(pointIndex[i]), new Integer(0));
			pointPowerDPvalue.put(Utils.makeKey(pointIndex[i]), new Integer(0));
			pointPowerEmpiricalAPvalue.put(Utils.makeKey(pointIndex[i]),
					new Integer(0));
			pointPowerEmpiricalDPvalue.put(Utils.makeKey(pointIndex[i]),
					new Integer(0));
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
		for (Iterator e = QTL.iterator(); e.hasNext();) {
			AbstractLoci al = (AbstractLoci) e.next();
			int chr = al.getChr()[0];
			double loc = al.getLocation()[0];
			ArrayList chromosome = (ArrayList) chromosomes.get(chr);
			for (int j = 0; j < chromosome.size(); j++) {
				double loc_map = ((Double) chromosome.get(j)).doubleValue();
				if (loc < loc_map) {
					chromosome.add(j, new Double(loc));
					al.setLoci(j);
					break;
				}
			}
		}
		distance = new double[chromosomes.size()][];
		for (int i = 0; i < chromosomes.size(); i++) {
			ArrayList chromosome = (ArrayList) chromosomes.get(i);
			distance[i] = new double[chromosome.size()];
			for (int j = 0; j < chromosome.size(); j++) {
				distance[i][j] = ((Double) chromosome.get(j)).doubleValue();
			}
		}
	}

	public abstract void correctLogBon();

	protected abstract void selectMarker() ;

	public void Simulation(Parameter1 p, ArrayList qtl, double[] env,
			double[][] w, double[][] m, ArrayList sm, boolean isSelectMarker) {
		Param1 = p;
		QTL = qtl;
		map = new double[m.length][];
		for (int i = 0; i < map.length; i++) {
			map[i] = new double[m[i].length];
			System.arraycopy(m[i], 0, map[i], 0, map[i].length);
		}
		weight = new double[w.length][];
		for (int i = 0; i < w.length; i++) {
			weight[i] = new double[w[i].length];
			System.arraycopy(w[i], 0, weight[i], 0, weight[i].length);
		}
		for (int i = 0; i < map.length; i++) {
			Num_interval += map[i].length == 1? map[i].length : map[i].length - 1;
		}
		if (sm != null) {
			selectedMarker = sm;
		} else {
			if (isSelectMarker) {
				selectMarker();
			}
		}
		correctLogBon();
		makeFullMap(QTL, map);
		setup(QTL, map, Param1);
		SimulationResults = new ArrayList();
		for (int i_rep = 0; i_rep < Param1.rep; i_rep++) {
			if (Param1.pt.compareTo("F2") == 0) {
				ap = new F2Population(Param1.populationSize,
						Param1.pheNum.length, Param1.pt, distance, Param1.seed
								+ i_rep, Param1.mu, env, Param1.sd, QTL,
						Param1.mf);
			} else if (Param1.pt.compareTo("DH") == 0) {
				ap = new DHPopulation(Param1.populationSize,
						Param1.pheNum.length, Param1.pt, distance, Param1.seed
								+ i_rep, Param1.mu, env, Param1.sd, QTL,
						Param1.mf);
			} else {
				ap = new BackCrossPopulation(Param1.populationSize,
						Param1.pheNum.length, Param1.pt, distance, Param1.seed
								+ i_rep, Param1.mu, env, Param1.sd, QTL,
						Param1.mf);
			}
			ap.ProducePopulation();

			for (int i = 0; i < Param1.pheNum.length; i++) {
				ap.ProducePhenotype(i, Param1.MU, Param1.T);
			}
			PrintStream cim = null;
			try {
				cim = new PrintStream(new BufferedOutputStream(
					new FileOutputStream("CIM.mcd")));
			} catch (Exception e) {
				
			}
			ap.CIMFormat(cim);
			cim.close();
			gs = new GenomeScan(ap, Param1.step);
			gs.CalculateIPP();

			if (Param1.permutation > 0 && i_rep == 0) {
				thresholdLOD = new ArrayList();
				thresholdAP = new ArrayList();
				thresholdDP = new ArrayList();
				thresholdPwald = new ArrayList();
				thresholdPF = new ArrayList();
				for (int i_permu = 0; i_permu < Param1.permutation; i_permu++) {
					ArrayList IMStatistic = MappingProcedure(i_rep, i_permu + 1);
					ArrayList IMLOD = new ArrayList();
					ArrayList IMAT = new ArrayList();
					ArrayList IMDT = new ArrayList();
					ArrayList IMAP = new ArrayList();
					ArrayList IMDP = new ArrayList();
					ArrayList IMPWald = new ArrayList();
					ArrayList IMPF = new ArrayList();
					int c = 0;
					int ct = 0;
					int bt = 0;

					for (int j = 0; j < IMStatistic.size(); j++) {
						PointMappingStatistic pms = (PointMappingStatistic) IMStatistic
								.get(j);
						IMLOD.add(new Double(pms.get_LOD()));
						IMAT.add(new Double(pms.get_tStatistic_additive()));
						IMDT.add(new Double(pms.get_tStatistic_dominance()));
						IMAP.add(new Double(pms.get_logP_additive()));
						IMPWald.add(new Double(pms.get_logP_wald()));
						IMPF.add(new Double(pms.get_logP_F()));
						if (pms.get_logP_additive() > LogBon) {
							c++;
						}
						if (pms.get_P_additive() < BonT) {
							bt++;
						}
						if (weight.length > 1) {
							IMDP.add(new Double(pms.get_logP_dominance()));
						}
					}
					if (c > 0) {
						typeI_LogBon += 1;
					}
					thresholdLOD.add(Collections.max(IMLOD));
					thresholdAP.add(Collections.max(IMAP));
					if (weight.length > 1) {
						thresholdDP.add(Collections.max(IMDP));
					}
					thresholdPF.add(Collections.max(IMPF));
					thresholdPwald.add(Collections.max(IMPWald));
				}
			}
			ArrayList IMLOD = MappingProcedure(i_rep, 0);
			SimulationResults.add(IMLOD);
		}
	}

	public abstract ArrayList MappingProcedure(int i_rep, int isPermutation);

	private void calculatePower() {
		ChiSquaredDistribution chi = new ChiSquaredDistributionImpl(weight.length);
		Collections.sort(thresholdAP);
		Collections.sort(thresholdDP);
		Collections.sort(thresholdLOD);
		Collections.sort(thresholdPwald);
		Collections.sort(thresholdPF);
		EmpiricalLogBonAP = ((Double) thresholdAP.get((new Double(
				(thresholdAP.size() - 1) * 0.95)).intValue()))
				.doubleValue();
		if(weight.length > 1) {
			EmpiricalLogBonDP = ((Double) thresholdDP.get((new Double(
					(thresholdDP.size() - 1) * 0.95)).intValue()))
					.doubleValue();
		}
		
		threshold_LOD = ((Double) thresholdLOD.get((new Double(
				(thresholdLOD.size() - 1) * 0.95)).intValue())).doubleValue();
		threshold_lodp = 0;
		double ln = threshold_LOD/0.217;
		threshold_Pwald = ((Double) thresholdPwald.get((new Double(
				(thresholdPwald.size() - 1) * 0.95)).intValue())).doubleValue();
		threshold_PF = ((Double) thresholdPF.get((new Double(
				(thresholdPF.size() - 1) * 0.95)).intValue())).doubleValue();
		try {
			threshold_lodp = -1*Math.log10(1-chi.cumulativeProbability(ln));
		} catch (Exception E) {
			E.printStackTrace(System.err);
		}
		System.out.println("BonT: " + BonT);
		System.out.println("LogBon: " + LogBon );
		System.out.println("EmpLogBonAP: " + EmpiricalLogBonAP);
		System.out.println("EmpLogBonDP: " + EmpiricalLogBonDP);
		System.out.println("threshold_LOD: " + "log10:" + threshold_LOD + " LR: " + ln + " P(ln) " + threshold_lodp );
		System.out.println("threshold_PWald: " + threshold_Pwald );
		System.out.println("threshold_PF: " + threshold_PF);
		for (Iterator e = SimulationResults.iterator(); e.hasNext();) {
			ArrayList result = (ArrayList) e.next();
			for (Iterator e1 = result.iterator(); e1.hasNext();) {
				PointMappingStatistic pms = (PointMappingStatistic) e1.next();
				String key = pms.getKey();
				if (Param1.permutation > 0 && pointPowerF.containsKey(key)
						&& pms.get_logP_F() > threshold_PF) {
					int p = ((Integer) pointPowerF.get(key)).intValue();
					p++;
					pointPowerF.put(key, new Integer(p));
				}			
				if (Param1.permutation > 0 && pointPowerWald.containsKey(key)
						&& pms.get_logP_wald() > threshold_Pwald) {
					int p = ((Integer) pointPowerWald.get(key)).intValue();
					p++;
					pointPowerWald.put(key, new Integer(p));
				}
				if (Param1.permutation > 0 && pointPowerLOD.containsKey(key)
						&& pms.get_LOD() > threshold_LOD) {
					int p = ((Integer) pointPowerLOD.get(key)).intValue();
					p++;
					pointPowerLOD.put(key, new Integer(p));
				}
				if (pointPowerAPvalue.containsKey(key)
						&& pms.get_logP_additive() > LogBon) {
					int p = ((Integer) pointPowerAPvalue.get(key)).intValue();
					p++;
					pointPowerAPvalue.put(key, new Integer(p));
				}
				if (pointPowerEmpiricalAPvalue.containsKey(key)
						&& pms.get_logP_additive() > EmpiricalLogBonAP) {
					int p = ((Integer) pointPowerEmpiricalAPvalue.get(key))
							.intValue();
					p++;
					pointPowerEmpiricalAPvalue.put(key, new Integer(p));
				}
				if (weight.length > 1) {
					if (pointPowerDPvalue.containsKey(key)
							&& pms.get_logP_dominance() > LogBon) {
						int p = ((Integer) pointPowerDPvalue.get(key))
								.intValue();
						p++;
						pointPowerDPvalue.put(key, new Integer(p));
					}
					if (pointPowerEmpiricalDPvalue.containsKey(key)
							&& pms.get_logP_dominance() > EmpiricalLogBonDP) {
						int p = ((Integer) pointPowerEmpiricalDPvalue.get(key))
								.intValue();
						p++;
						pointPowerEmpiricalDPvalue.put(key, new Integer(p));
					}
				}
			}
		}
		Set keys = pointPowerLOD.keySet();
		System.out.println("type I LogBon: " + typeI_LogBon
				/ Param1.permutation);
		System.out.println("QTL\tFP\tWald\tLOD\tBonAP\tEmpAP\tBonDP\tEmpDP");
		for (Iterator e = keys.iterator(); e.hasNext();) {
			String key = (String) e.next();
			System.out.print(key + "\t");
			System.out.print((Integer) pointPowerF.get(key) + "\t");
			System.out.print((Integer) pointPowerWald.get(key) + "\t");
			System.out.print((Integer) pointPowerLOD.get(key) + "\t");
			System.out.print((Integer) pointPowerAPvalue.get(key) + "\t");
			System.out.print((Integer) pointPowerEmpiricalAPvalue.get(key)+ "\t");
			if (weight.length > 1) {
				System.out.print(pointPowerDPvalue.get(key) + "\t");
				System.out.println(pointPowerEmpiricalDPvalue.get(key) + "\t");
			} else {
				System.out.println("--\t--\t");
			}
		}
	}

	public double get_threshold_LOD() {
		return threshold_LOD;
	}
	
	public double get_threshold_LODP() {
		return threshold_lodp;
	}
	
	public double get_threshold_PWald() {
		return threshold_Pwald;
	}
	
	public double get_threshold_PF() {
		return threshold_PF;
	}

	public void SummuarySimulation() {
//		calculatePower();
		PrintStream Pout1 = null;
		PrintStream Pout2 = null;
		PrintStream Pout3 = null;
		PrintStream Pout4 = null;
		PrintStream Pout5 = null;
		PrintStream Pout6 = null;
		ChiSquaredDistribution chi = new ChiSquaredDistributionImpl(weight.length);
		try {
			Pout1 = new PrintStream(new BufferedOutputStream(
					new FileOutputStream("LogAP.txt")));
			Pout2 = new PrintStream(new BufferedOutputStream(
					new FileOutputStream("LNP.txt")));
			Pout3 = new PrintStream(new BufferedOutputStream(
					new FileOutputStream("LN.txt")));
			Pout4 = new PrintStream(new BufferedOutputStream(
					new FileOutputStream("Wald.txt")));
			Pout5 = new PrintStream(new BufferedOutputStream(
					new FileOutputStream("F.txt")));
			Pout6 = new PrintStream(new BufferedOutputStream(
					new FileOutputStream("LogDP.txt")));
		} catch (Exception E) {
			E.printStackTrace(System.err);
		}
		for (int i = 0; i < SimulationResults.size(); i++) {
			ArrayList PointStatistics = (ArrayList) SimulationResults.get(i);
			double[] LRT = new double[PointStatistics.size()];
			double[] WALD = new double[PointStatistics.size()]; 
			DescriptiveStatistics dsLRT = new DescriptiveStatisticsImpl();
			DescriptiveStatistics dsWALD = new DescriptiveStatisticsImpl();
			SimpleRegression sr = new SimpleRegression();
			for (int j = 0; j < PointStatistics.size(); j++) {
				PointMappingStatistic pms = (PointMappingStatistic) PointStatistics
						.get(j);
				double lod = pms.get_LOD()>0 ? pms.get_LOD() : 0;
				double ln = pms.get_LN();
				double p = 0;
				try {
					p = chi.cumulativeProbability(ln);
				} catch (Exception E) {
					E.printStackTrace(System.err);
				}
				double logLOD = -1 * Math.log10(1-p);
				Pout1.print(pms.get_logP_additive() + " ");
				Pout2.print(logLOD + " ");
				Pout3.print(pms.get_LN() + " ");
				Pout4.print(pms.get_wald() + " ");
				Pout5.print(pms.get_logP_F() + " ");
				Pout6.print(pms.get_logP_dominance() + " ");

				dsLRT.addValue(pms.get_LN());
				dsWALD.addValue(pms.get_wald());
				LRT[j] = pms.get_LN();
				WALD[j] = pms.get_wald();
				sr.addData(pms.get_LN(), pms.get_wald());
			}
			System.out.println(dsLRT.getMean() + " +- " + dsLRT.getStandardDeviation() + " " + dsWALD.getMean() + " +- " + dsWALD.getStandardDeviation() + " cor " + sr.getR());
			dsLRT.clear();
			dsWALD.clear();
			sr.clear();
			Pout1.println();
			Pout2.println();
			Pout3.println();
			Pout4.println();
			Pout5.println();
			Pout6.println();
		}
		Pout1.close();
		Pout2.close();
		Pout3.close();
		Pout4.close();
		Pout5.close();
		Pout6.close();
	}
}
