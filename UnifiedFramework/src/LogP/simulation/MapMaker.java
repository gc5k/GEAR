package LogP.simulation;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;

public class MapMaker {
	private int chr;
	private int[] marker;
	private double[][] map;
	private Random rnd;
	public MapMaker() {
		rnd = new Random();
	}

	public void MakeRandomMap(double[] length, int[] markers, long seed, double[] criterion) {
		rnd.setSeed(seed);
		map = new double[length.length][];
		for (int i = 0; i < length.length; i++) {
			Set m = new TreeSet();
			m.add(new Double(0));
			m.add(new Double(length[i]));
			int idx = 2;
			while(idx < markers[i]) {
				double point = rnd.nextFloat() * length[i];
				double min = 1;
				for (Iterator e = m.iterator(); e.hasNext(); ) {
					double p = ((Double) e.next()).doubleValue();
					if (Math.abs(p - point) < min ) {
						min = Math.abs(p-point);
					}
				}
				if (min > criterion[0] && min < criterion[1]) {
					m.add(new Double (point));
					idx++;
				}
			}
			map[i] = new double[markers[i]];
			int c = 0;
			for (Iterator e = m.iterator(); e.hasNext(); ) {
				map[i][c++] = ((Double) e.next());
			}
		}
	}

	public double[][] getMap() {
		return map;
	}

	public ArrayList getRandomMarker(double[][] m, long seed) {
		rnd.setSeed(seed);
		ArrayList sm = new ArrayList(); 
		for (int i = 0; i < m.length; i++) {
			int c = rnd.nextInt(m[i].length-1)+1;
			Set rm = new TreeSet();
			while(rm.size() < c) {
				boolean flag = true;
				int mi = rnd.nextInt(m[i].length);
				for (Iterator e = rm.iterator(); e.hasNext();) {
					int mis = ((Integer) (e.next())).intValue();
					if (mis == mi) {
						flag = false;
						break;
					}
				}
				if(flag) {
					rm.add(new Integer(mi));
				}
			}
			for (Iterator e = rm.iterator(); e.hasNext(); ) {
				ArrayList point = new ArrayList();
				point.add(new Integer(i));
				point.add(e.next());
				sm.add(point);
			}
		}
		return sm;
	}

	public ArrayList getRandomMarker(long seed) {
		rnd.setSeed(seed);
		ArrayList sm = new ArrayList(); 
		for (int i = 0; i < map.length; i++) {
			int c = rnd.nextInt(map[i].length-1)+1;
			Set rm = new TreeSet();
			while(rm.size() < c) {
				boolean flag = true;
				int mi = rnd.nextInt(map[i].length);
				for (Iterator e = rm.iterator(); e.hasNext();) {
					int mis = ((Integer) (e.next())).intValue();
					if (mis == mi) {
						flag = false;
						break;
					}
				}
				if(flag) {
					rm.add(new Integer(mi));
				}
			}
			for (Iterator e = rm.iterator(); e.hasNext(); ) {
				ArrayList point = new ArrayList();
				point.add(new Integer(i));
				point.add(e.next());
				sm.add(point);
			}
		}
		return sm;
	}

	public static void main(String[] args) throws IOException {
		MapMaker mm = new MapMaker();
		double[] len = {2.0};
		int[] marker = {20};
		long seed = 1;
		double[] criterion = {0.05, 0.2};
		mm.MakeRandomMap(len, marker, seed, criterion);
		double[][] map = mm.getMap();
		for (int i = 0; i < map.length; i++) {
			for (int j = 0; j < map[i].length; j++) {
				System.out.print(new java.text.DecimalFormat("0.00").format(map[i][j]) + " ");
			}
			System.out.println();
		}
		for (int i = 200; i < 230; i++) {
			ArrayList sm = mm.getRandomMarker(i);
			for (Iterator e = sm.iterator(); e.hasNext(); ) {
				ArrayList point = (ArrayList) e.next();
				System.out.print(point.get(0) + "->" + point.get(1) + ",");
			}
			System.out.println();
		}
	}
}
