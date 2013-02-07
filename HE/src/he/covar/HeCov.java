package he.covar;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;
import java.util.TreeSet;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;

import parameter.Parameter;

import util.FileProcessor;
import util.NewIt;

public class HeCov {
	private final String delim = "\\s+";
	private double[][] y;
	private boolean[] cov_flag;

	private boolean[] flag;

	private String qcov_file;
	private int[] qcov_idx;
	private double[][] q_cov;

	private String cov_file;
	private int[] cov_idx;
	private ArrayList<ArrayList<String>> cov;
	private ArrayList<HashMap<String, Integer>> cov_class;

	private HashMap<String, Integer> ID2Idx;

	public HeCov(double[][] p, boolean[] f, HashMap<String, Integer> I, String qf, int[] q_i, String cf, int[] c_i) {
		
		y = p;

		flag = f;
		ID2Idx = I;
		qcov_file = qf;
		qcov_idx = q_i;

		cov_file = cf;
		cov_idx = c_i;

		cov_flag = new boolean[flag.length];
		Arrays.fill(cov_flag, false);
		generate_Res();
	}

	private void generate_Res() {

		if( qcov_file != null ) {
			
			BufferedReader reader = FileProcessor.FileOpen(qcov_file);
			String line;
			ArrayList<String> c = null;
			int _c = 0;
			try {
				while ((line = reader.readLine()) != null) {

					String[] s = line.split(delim);
					StringBuilder sb = new StringBuilder();
					sb.append(s[0] + "." + s[1]);

					if (_c==0) {
						if(qcov_idx == null) {//use all the covariates if no covariates are selected
							qcov_idx = new int[s.length-2];
							for(int i = 0; i < qcov_idx.length; i++) {
								qcov_idx[i] = i+1;
							}
						} else {
							for(int i = 0; i < qcov_idx.length; i++) {
								System.out.println("qcov" + qcov_idx[i]);
							}
						}
						q_cov = new double[flag.length][qcov_idx.length];
					}

					if (ID2Idx.containsKey(sb.toString())) {
						boolean f = true;
						int idx = ID2Idx.get(sb.toString()).intValue();
						c = NewIt.newArrayList();
						for (int i = 0; i < qcov_idx.length; i++) {
							if(Parameter.isNA(s[1+qcov_idx[i]])) {
								f = false;
							}
							c.add(s[1+qcov_idx[i]]);
						}
						
						cov_flag[idx] = f;
						if(f) {
							for (int i = 0; i < c.size(); i++) {
								q_cov[idx][i] = Double.parseDouble(c.get(i));
							}
						}
					}
					_c++;
				}
				reader.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		if( cov_file != null ) {

			BufferedReader reader = FileProcessor.FileOpen(cov_file);
			String line;
			ArrayList<String> c = null;
			int _c = 0;
			try {
				while ((line = reader.readLine()) != null) {
					String[] s = line.split(delim);
					StringBuilder sb = new StringBuilder();
					sb.append(s[0] + "." + s[1]);

					if (_c==0) {
						if (cov_idx == null) {//use all the covariates
							cov_idx = new int[s.length-2];
							for (int i = 0; i < cov_idx.length; i++) {
								cov_idx[i] = i+1;
							}
						}

						cov = NewIt.newArrayList();
						for (int i = 0; i < flag.length; i++) {
							ArrayList<String> S = NewIt.newArrayList();
							S.ensureCapacity(cov_idx.length);
							cov.add(S);
						}
						cov_class = NewIt.newArrayList();
						for (int i = 0; i < cov_idx.length; i++) {
							HashMap<String, Integer> hm = NewIt.newHashMap();
							cov_class.add(hm);
						}
					}

					if (ID2Idx.containsKey(sb.toString())) {
						boolean f = true;
						int idx = ID2Idx.get(sb.toString()).intValue();
						c = NewIt.newArrayList();
						for (int i = 0; i < cov_idx.length; i++) {
							if(Parameter.isNA(s[1+cov_idx[i]])) {
								f = false;
							} else {// exclude NA individuals
								c.add(s[1+cov_idx[i]]);
							}
						}
						
						if ( qcov_file != null ) {//should be existed in both
							cov_flag[idx] &= f;
						} else {//if exists in this file
							cov_flag[idx] = f;				
						}
						if (cov_flag[idx]) {
							cov.set(idx, c);
							for (int i = 0; i < c.size(); i++) {
								String s_ = s[1+cov_idx[i]];
								HashMap<String, Integer> hm = cov_class.get(i);
								if (!hm.containsKey(s_)) {
									int size = hm.size() + 1;
									hm.put(s_, size);
								}
							}
						}
					}
					_c++;
				}
				reader.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		int dim = 0;
		for (int i = 0; i < cov_flag.length; i++) {
			if(cov_flag[i]) {
				dim++;
			} else {
				System.out.println(i);
			}
		}

		System.out.println("dimension: " + dim);
		double[][] phe = new double[dim][1];

		int cn = 0;
		for (int i = 0; i < cov_flag.length; i++) {
			if (!cov_flag[i]) continue;
			phe[cn][0] = y[i][1];
			cn++;
		}

		double[][] x1 = null;
		if (qcov_file != null) {
			x1 = new double[dim][qcov_idx.length];

			int cn1 = 0;
			for (int i = 0; i < cov_flag.length; i++) {
				if (!cov_flag[i]) continue;
				for (int j = 0; j < q_cov[i].length; j++) {
					x1[cn1][j] = q_cov[i][j];					
				}
				cn1++;
			}
		}
		
		double[][] x2 = null;
		if(cov_file != null) {
			x2 = new double[dim][cov_idx.length];

			int cn2 = 0;
			for (int i = 0; i < cov_flag.length; i++) {
				if (!cov_flag[i]) continue;
				ArrayList<String> c = cov.get(i);
				for (int j = 0; j < c.size(); j++) {
					HashMap<String, Integer> hm = cov_class.get(j);
					int code = hm.get(c.get(j)).intValue();
					x2[cn2][j] = code;
				}
				cn2++;
			}
		}

		double[][] x3 = null;
		if(qcov_file!=null && cov_file != null) {
			x3 = new double[cn][qcov_idx.length+cov_idx.length+1];
			for(int i = 0; i < dim; i++) {
				x3[i][0] = 1;
				System.arraycopy(x1[i], 0, x3[i], 1, qcov_idx.length);
				System.arraycopy(x2[i], 0, x3[i], 1+qcov_idx.length, cov_idx.length);
			}
		} else if(qcov_file != null) {
			x3 = new double[cn][qcov_idx.length+1];
			for(int i = 0; i < dim; i++) {
				x3[i][0] = 1;
				System.arraycopy(x1[i], 0, x3[i], 1, qcov_idx.length);
			}			
		} else if(cov_file != null) {
			x3 = new double[cn][cov_idx.length+1];
			for(int i = 0; i < dim; i++) {
				x3[i][0] = 1;
				System.arraycopy(x2[i], 0, x3[i], 1, cov_idx.length);
			}
		}

		System.out.println(x3.length + " " + x3[0].length);
		RealMatrix Y = new Array2DRowRealMatrix(phe);
		RealMatrix X = new Array2DRowRealMatrix(x3);

		RealMatrix XtX = X.transpose().multiply(X);
		RealMatrix XtY = X.transpose().multiply(Y);

		RealMatrix XtX_Inv = (new LUDecompositionImpl(XtX)).getSolver().getInverse();
		RealMatrix Mat_B = XtX_Inv.multiply(XtY);

		System.out.println(Mat_B);
		RealMatrix Pre = X.multiply(Mat_B);

		int cc=0;
		for(int i = 0; i < cov_flag.length; i++) {
			flag[i] &= cov_flag[i];
			if (!cov_flag[i]) continue;
			y[i][1] -= Pre.getEntry(cc++, 0);
		}
	}
}
