package simulation;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;

import org.apache.commons.math.MathException;
import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.random.RandomDataImpl;

import parameter.Parameter;

public class SimuFamily {
	private RandomDataImpl rnd;
	private long seed = 2011;

	private final String[] A={"A", "C"};
	private int NFam = 100;
	private int NMarker = 10;
	private int[] NKid = null;
	private int[] NAffKid = null;
	private double[] LD = null;
	private double[] rec = null;
	private double[] maf = null;
	private double[] effect = null;
	
	private PrintWriter pedout = null;
	private PrintWriter map = null;
	private PrintWriter phe = null;
	private PrintWriter cov = null;
	
	private String out = "simuFam";
	private Parameter par;

	public SimuFamily(Parameter p) {

		par = p;
		NFam = par.simu_fam_size;
		NMarker = par.simu_fam_marker;
		seed = par.seed;

		initial();
	}

	public SimuFamily() {

		initial();

	}

	private void initial() {

		rnd = new RandomDataImpl();
		rnd.reSeed(seed);

		NKid = new int[NFam];
		Arrays.fill(NKid, 2);
		NAffKid = new int[NFam];
		Arrays.fill(NAffKid, 1);

		maf = new double[NMarker];
		Arrays.fill(maf, 0.5);
		LD = new double[NMarker];
		Arrays.fill(LD, 0);
		rec = new double[NMarker];
		Arrays.fill(rec, 0.5);
		rec[0] = maf[0];

	}
	
	public void generateSample() {

		try {
			pedout = new PrintWriter(new BufferedWriter(new FileWriter(out
					+ ".ped")));
			map = new PrintWriter(new BufferedWriter(new FileWriter(out
					+ ".map")));
			phe = new PrintWriter(new BufferedWriter(new FileWriter(out
					+ ".phe")));
			cov = new PrintWriter(new BufferedWriter(new FileWriter(out
					+ ".cov")));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		for (int i = 0; i < NFam; i++) {
			generateNuclearFamily(NKid[i], NAffKid[i], i+1);
		}
		
		writeMap();
		
		pedout.close();
		map.close();
		phe.close();
		cov.close();

	}

	private void generateNuclearFamily(int nkid, int affKid, int famIdx) {
		
		RealMatrix pchr = sampleChromosome();
		RealMatrix mchr = sampleChromosome();
		RealMatrix[] kchr = new RealMatrix[nkid];
		for(int i = 0; i< nkid; i++ ) {
			kchr[i] = generateBaby(pchr, mchr);
		}

		writeFile(pchr, mchr, kchr, famIdx);

	}

	private RealMatrix sampleChromosome() {

		double[] g = new double[maf.length];
		double[][] v = new double[maf.length][2];
		for (int i = 0; i < maf.length; i++) {
			for (int j = 0; j < 2; j++) {
				double r = rnd.nextUniform(0, 1);
				if (i == 0) {
					v[i][j] = r < maf[i] ? 0 : 1;
				} else {
					double d = rnd.nextUniform(0, 1);
					int a = (int) v[i - 1][j];
					double f1 = a == 0 ? maf[i - 1] : (1 - maf[i - 1]);
					double f2 = a == 0 ? maf[i] : (1 - maf[i]);
					v[i][j] = d < (f1 * f2 + LD[i - 1]) / f1 ? v[i - 1][j]
							: (1 - v[i - 1][j]);
				}
			}
			g[i] = v[i][0] + v[i][1];
		}

		RealMatrix chr = new Array2DRowRealMatrix(v);
		return chr;
	}

	private RealMatrix generateBaby(RealMatrix pchr, RealMatrix mchr) {
		double[][] v = new double[maf.length][2];

		for (int i = 0; i < 2; i++) {
			RealMatrix chr = i == 0 ? pchr:mchr;
			int idx = 1;
			for (int j = 0; j < maf.length; j++) {
				double r = rnd.nextUniform(0, 1);
				idx = r < rec[j] ? 1 - idx : idx;
				v[j][i] = chr.getEntry(j, idx);
			}
		}

		RealMatrix chr = new Array2DRowRealMatrix(v);
		return chr;
	}

	public void writeFile(RealMatrix pchr, RealMatrix mchr, RealMatrix[] kchr, int famIdx) {

		int fid=famIdx * 10000;
		int pid = fid+1;
		int mid = fid+2;

		pedout.print(fid + " ");
		pedout.print(pid + " ");
		pedout.print(0 + " ");
		pedout.print(0 + " ");
		pedout.print(1 + " ");
		pedout.print(1 + " ");

		for (int i = 0; i < pchr.getRowDimension(); i++) {
			pedout.print(A[(int) (pchr.getEntry(i, 0))] + " " +  A[(int) pchr.getEntry(i, 1)] + " ");
		}
		pedout.print("\n");

		pedout.print(fid + " ");
		pedout.print(mid + " ");
		pedout.print(0 + " ");
		pedout.print(0 + " ");
		pedout.print(2 + " ");
		pedout.print(1 + " ");

		for (int i = 0; i < mchr.getRowDimension(); i++) {
			pedout.print(A[(int) (mchr.getEntry(i, 0))] + " " +  A[(int) mchr.getEntry(i, 1)] + " ");
		}
		pedout.print("\n");

		for (int i = 0; i < kchr.length; i++) {
			pedout.print(fid + " ");
			pedout.print((fid+3+i) + " ");
			pedout.print(pid + " ");
			pedout.print(mid + " ");
			pedout.print(2 + " ");
			try {
				pedout.print((rnd.nextBinomial(1, 0.5) + 1)+ " ");
			} catch (MathException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

			for (int j = 0; j < kchr[i].getRowDimension(); j++) {
				pedout.print(A[(int) (kchr[i].getEntry(j, 0))] + " " +  A[(int) kchr[i].getEntry(j, 1)] + " ");
			}
			pedout.print("\n");
		}
	}

	private void writeMap() {

		for (int i = 0; i < maf.length; i++) {
			map.print(1 + " ");
			map.print("rs" + i + " ");
			map.print(i / (maf.length * 1.0) + " ");
			map.println(i * 100);
		}

	}

	public static void main(String[] args) {
		SimuFamily sf = new SimuFamily();
		sf.generateSample();
	}
}
