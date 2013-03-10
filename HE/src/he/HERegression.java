package he;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;
import java.util.zip.GZIPInputStream;

import parameter.Parameter;
import util.FileProcessor;

import org.apache.commons.math.MathException;
import org.apache.commons.math.linear.*;
import org.apache.commons.math.distribution.*;

public class HERegression {

	private final String delim = "\\s+";
	private boolean[] flag;
	HashMap<String, Integer> ID;
	private String grmFile;
	private String grmID;
	private String keepFile;
	private String phenoFile;
	private String output;
	private int[] mpheno;

	private boolean reverse;
	private boolean k_button;
	private double k;
	private boolean[] heType;

	private HashMap<String, Integer> ID2Idx;

	private double yyProd;

	private double[][] XtX;
	private double[] XtY;
	private int dim;
	private double P;
	private boolean isCC = false;

	StringBuffer sb = new StringBuffer();

	public HERegression() {
		grmFile = Parameter.INSTANCE.grm;
		grmID = Parameter.INSTANCE.grm_id;
		keepFile = Parameter.INSTANCE.keepFile;
		phenoFile = Parameter.INSTANCE.pheno;
		mpheno = Parameter.INSTANCE.mpheno;
		reverse = Parameter.INSTANCE.reverse;
		k_button = Parameter.INSTANCE.k_button;
		k = Parameter.INSTANCE.k;
		output = Parameter.INSTANCE.out;
		heType = Parameter.INSTANCE.heType;

		XtX = new double[mpheno.length + 1][mpheno.length + 1];
		XtY = new double[mpheno.length + 1];

		String line;

		// *********************************** read grm id
		BufferedReader reader = FileProcessor.FileOpen(grmID);

		int i2 = 0;
		ID2Idx = new HashMap<String, Integer>();
		try {
			while ((line = reader.readLine()) != null) {
				String[] s = line.split(delim);
				StringBuilder sb = new StringBuilder();
				sb.append(s[0] + "." + s[1]);
				ID2Idx.put(sb.toString(), i2++);
			}
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		flag = new boolean[ID2Idx.size()];
		Arrays.fill(flag, false);

		// *********************************** read pheno file
		reader = FileProcessor.FileOpen(phenoFile);

		double[][] y = new double[flag.length][mpheno.length + 1];

		try {
			while ((line = reader.readLine()) != null) {
				String[] s = line.split(delim);
				StringBuilder sb = new StringBuilder();
				sb.append(s[0] + "." + s[1]);
				if (ID2Idx.containsKey(sb.toString())) {
					int ii = ID2Idx.get(sb.toString());
					boolean f = true;
					y[ii][0] = 1;
					for (int j = 0; j < mpheno.length; j++) {
						if (Parameter.isNA(s[1 + mpheno[j]])) {
							f = false;
							break;
						} else {
							y[ii][j + 1] = Double.parseDouble(s[1 + mpheno[j]]);
						}
					}
					flag[ii] = f;
				}
			}
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		// ************************keep
		if (keepFile != null) {
			File keepF = new File(keepFile);
			reader = FileProcessor.FileOpen(keepFile);
			boolean[] ff = new boolean[flag.length];
			Arrays.fill(ff, false);
			try {
				while ((line = reader.readLine()) != null) {
					String[] s = line.split(delim);
					StringBuilder sb = new StringBuilder(s[0] + "." + s[1]);
					if (ID2Idx.containsKey(sb.toString())) {
						int ii = ID2Idx.get(sb.toString());
						ff[ii] = true;
					}
				}
				reader.close();
			} catch (IOException e) {
				e.printStackTrace();
			}

			for (int i = 0; i < ff.length; i++) {
				flag[i] &= ff[i];
			}
		}

		int Len = 0;
		for (int i = 0; i < flag.length; i++)
			if (flag[i]) Len++;
		dim = Len*(Len-1)/2;
		
		// ************************************standardising
		double[] ss = new double[y[0].length-1];
		double[] ssx = new double[y[0].length-1];
		for (int i = 0; i < flag.length; i++) {
			if(!flag[i]) continue; 
			for(int j = 0; j < ss.length; j++) {
				ss[j] += y[i][j+1];
				ssx[j] += y[i][j+1] * y[i][j+1];
			}
		}
		double[] sd = new double[ssx.length];
		for (int i = 0; i < sd.length; i++) {
			ss[i] /= Len;
			sd[i] = Math.sqrt((ssx[i] - Len * ss[i] * ss[i])/(Len-1));
		}

		if (Parameter.INSTANCE.scale) {
			System.out.println("standardising phentoype.");
			for (int i = 0; i < flag.length; i++) {
				if(!flag[i]) continue;
				for (int j = 1; j < y[i].length; j++) {
					y[i][j] = (y[i][j] - ss[j-1])/sd[j-1];
				}
			}
		}

		// *************************************read grm file
		FileInputStream fin = null;
		try {
			fin = new FileInputStream(grmFile);
		} catch (FileNotFoundException e1) {
			e1.printStackTrace();
		}
		GZIPInputStream gzis = null;
		try {
			gzis = new GZIPInputStream(fin);
		} catch (IOException e1) {
			e1.printStackTrace();
		}
		InputStreamReader xover = new InputStreamReader(gzis);
		BufferedReader is = new BufferedReader(xover);

		if (reverse) {
			try {
				while ((line = is.readLine()) != null) {
					String[] s = line.split(delim);
					int id1 = Integer.parseInt(s[0]) - 1;
					int id2 = Integer.parseInt(s[1]) - 1;
					if (id1 == id2)
						continue;
					if (!(flag[id1] & flag[id2]))
						continue;
					for (int j = 0; j < y[id1].length; j++) {
						double Xj = (y[id1][j] - y[id2][j])
								* (y[id1][j] - y[id2][j]);
						for (int k = j; k < y[id2].length; k++) {
							if (k == 0 && j == 0)
								XtX[j][k] += 1;
							double Xk = (y[id1][k] - y[id2][k])
									* (y[id1][k] - y[id2][k]);
							if (j == 0) {
								XtX[j][k] += Xk;
							} else {
								XtX[j][k] += Xj * Xk;
							}
						}
						if (j == 0) {
							XtY[j] += Double.parseDouble(s[3]);
						} else {
							XtY[j] += Xj * Double.parseDouble(s[3]);
						}
					}
					yyProd += Double.parseDouble(s[3])
							* Double.parseDouble(s[3]);
				}
				is.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		} else {
			XtX = new double[2][2];
			XtY = new double[2];
			try {
				HashMap<Double, Integer> cat = new HashMap<Double, Integer>();
				while ((line = is.readLine()) != null) {
					String[] s = line.split(delim);
					int id1 = Integer.parseInt(s[0]) - 1;
					int id2 = Integer.parseInt(s[1]) - 1;
					if (id1 == id2)
						continue;
					if (!(flag[id1] & flag[id2]))
						continue;
					double ds = 0;
					if (heType[Parameter.he_sd]) {
						ds = (y[id1][1] - y[id2][1]) * (y[id1][1] - y[id2][1]);
					} else if(heType[Parameter.he_ss]) {
						ds = (y[id1][1] + y[id2][1]) * (y[id1][1] + y[id2][1]);
					} else if(heType[Parameter.he_cp]) {
						ds = y[id1][1] * y[id2][1];
					}

					if (cat.containsKey(y[id1][1])) {
						Integer I = (Integer) cat.get(y[id1][1]);
						I++;
						cat.put(y[id1][1], I);
					} else {
						cat.put(y[id1][1], 1);
					}
					if (cat.containsKey(y[id2][1])) {
						Integer I = (Integer) cat.get(y[id2][1]);
						I++;
						cat.put(y[id2][1], I);
					} else {
						cat.put(y[id2][1], 1);
					}
					yyProd += ds * ds;
					XtX[0][0]++;
					double g = Double.parseDouble(s[3]);
					XtX[0][1] += g;
					XtX[1][0] += g;
					XtX[1][1] += g * g;

					XtY[0] += ds;
					XtY[1] += g * ds;
				}

				if (cat.size() == 2) {
					isCC = true;
					Set<Double> set = cat.keySet();
					Iterator<Double> it = set.iterator();
					Double k1 = it.next();
					Double k2 = it.next();
					Integer c1 = cat.get(k1);
					Integer c2 = cat.get(k2);
					
					if(k1 > k2) {
						P = c1.doubleValue()/(c1.doubleValue() + c2.doubleValue());
					} else {
						P = c2.doubleValue()/(c1.doubleValue() + c2.doubleValue());
					}
				}
				is.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		for (int i = 0; i < XtX[0].length; i++) {
			for (int j = 0; j < i; j++) {
				XtX[i][j] = XtX[j][i];
			}
		}
	}

	public void Regression() {
		DecimalFormat fmt = new DecimalFormat("#.###E0");

		RealMatrix Mat_XtX = new Array2DRowRealMatrix(XtX);
		RealMatrix Mat_XtY = new Array2DRowRealMatrix(XtY);

		RealMatrix Mat_XtX_Inv = (new LUDecompositionImpl(Mat_XtX)).getSolver().getInverse();
		RealMatrix Mat_B = Mat_XtX_Inv.multiply(Mat_XtY);
		RealMatrix Mat_K = Mat_B.transpose().multiply(Mat_XtY);

		double mse = (yyProd - Mat_K.getEntry(0, 0))
				/ (dim - mpheno.length - 1);
		RealMatrix v = Mat_XtX_Inv.scalarMultiply(mse);

		sb.append("HE mode: ");
		if (heType[Parameter.he_sd]) {
			sb.append("squared difference (yi-yj)^2\n");
		} else if(heType[Parameter.he_ss]) {
			sb.append("squared sum (yi+yj)^2\n");
		} else if(heType[Parameter.he_cp]) {
			sb.append("cross-product [yi-E(y)][yj-E(y)]\n");
		}
		sb.append("grm: " + grmFile + "\n");
		sb.append("grm id: " + grmID + "\n");
		sb.append("keep list: " + keepFile + "\n");
		sb.append("pheno: " + phenoFile + "\n");
		sb.append("mpheno: ");
		for (int i = 0; i < mpheno.length; i++) {
			sb.append(mpheno[i] + " ");
		}
		sb.append("\n");
		
		sb.append("reverse: " + reverse + "\n");
		if(k_button) {
			sb.append("k: " + k + "\n");
		}
		sb.append("Coef\t" + "Estimate \t" + "se" + "\n");
		
		for (int i = 0; i < Mat_B.getRowDimension(); i++) {
			sb.append("b" + i + "\t" + fmt.format(Mat_B.getEntry(i, 0)) + "\t" + fmt.format(Math.sqrt(v.getEntry(i, i))) + "\n");
		}

		if (!reverse) {
			double h_o = 0;
			if (heType[Parameter.he_sd]) {
				h_o = Mat_B.getEntry(1, 0) / Mat_B.getEntry(0, 0) * (-1);
			} else if (heType[Parameter.he_ss]) {
				h_o = Mat_B.getEntry(1, 0) / Mat_B.getEntry(0, 0);
			} else if (heType[Parameter.he_cp]) {
				h_o = Mat_B.getEntry(1, 0);
			}

			double u_b0 = Mat_B.getEntry(0, 0);
			double u_b1 = Mat_B.getEntry(1, 0);

			double v_b0 = v.getEntry(0, 0);
			double v_b1 = v.getEntry(1, 1);

			double v_ho = (u_b1/u_b0) * (u_b1/u_b0) * (v_b0/(u_b0 * u_b0) + v_b1/(u_b1 * u_b1) - 2*v.getEntry(0, 1)/(u_b0*u_b1));
			if (heType[Parameter.he_cp]) {
				v_ho = Math.sqrt(v_b1);
			}
			sb.append("h2(o): " + fmt.format(h_o) + "\t" + fmt.format(v_ho) + "\n");

			if (k_button && isCC) {
				NormalDistributionImpl Norm = new NormalDistributionImpl();
				double q = 0;
				try {
					q = Norm.inverseCumulativeProbability(1-k);
				} catch (MathException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				double z = 1 / (Math.sqrt(2*3.1416926)) * Math.exp(-q*q/2);
				double h_l = h_o * k * (1-k) * k * (1-k) / (z * z * P * (1-P));

				double f = (k * (1-k) * k * (1-k))/( z * z * P * (1-P));
				double v_hl = v_ho * f * f;
				sb.append("h2(l): " + fmt.format(h_l) + "\t" + fmt.format(v_hl) + "\n");
			}
		}

		sb.append("\n");
		sb.append("variance-covariance\n");
		for (int i = 0; i < v.getRowDimension(); i++) {
			for (int j = 0; j <= i; j++) {
				sb.append(fmt.format(v.getEntry(i, j)) + " ");
			}
			sb.append("\n");
		}
				
		System.out.println(sb);

		if(output != null) {
			StringBuilder fsb = new StringBuilder();
			fsb.append(output);
			fsb.append(".he");
			File of = new File(fsb.toString());
			PrintWriter pw = null;
			try {
				pw = new PrintWriter(of);
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			pw.append(sb);
			pw.close();
		}
		
	}

	public static void main(String[] args) {
		Parameter.INSTANCE.commandListener(args);
		HERegression HER = new HERegression();
		HER.Regression();
	}

}
