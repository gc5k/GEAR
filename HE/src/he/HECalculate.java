package he;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;
import java.util.zip.GZIPInputStream;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistributionImpl;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.Array2DRowRealMatrix;

import he.endian.LittleEndianDataInputStream;

import parameter.Parameter;

public class HECalculate {

	private final String delim = "\\s+";
	private HERead heReader;
	private RealMatrix Mat_B;

	public HECalculate(HERead h) {
		heReader = h;
		Calculate();
	}

	public void Calculate() {

		heReader.lambda = new Lambda();

		String line;
		if (heReader.reverse) {
			// try {
			// while ((line = heReader.is.readLine()) != null) {
			// String[] s = line.split(delim);
			// int id1 = Integer.parseInt(s[0]) - 1;
			// int id2 = Integer.parseInt(s[1]) - 1;
			// if (id1 == id2)
			// continue;
			// if (!(heReader.flag[id1] & heReader.flag[id2]))
			// continue;
			// for (int j = 0; j < heReader.y[id1].length; j++) {
			// double Xj = (heReader.y[id1][j] - heReader.y[id2][j])
			// * (heReader.y[id1][j] - heReader.y[id2][j]);
			// for (int k = j; k < heReader.y[id2].length; k++) {
			// if (k == 0 && j == 0)
			// heReader.XtX[j][k] += 1;
			// double Xk = (heReader.y[id1][k] - heReader.y[id2][k])
			// * (heReader.y[id1][k] - heReader.y[id2][k]);
			// if (j == 0) {
			// heReader.XtX[j][k] += Xk;
			// } else {
			// heReader.XtX[j][k] += Xj * Xk;
			// }
			// }
			// if (j == 0) {
			// heReader.XtY[j] += Double.parseDouble(s[3]);
			// } else {
			// heReader.XtY[j] += Xj * Double.parseDouble(s[3]);
			// }
			// }
			// heReader.yyProd += Double.parseDouble(s[3])
			// * Double.parseDouble(s[3]);
			// }
			// heReader.is.close();
			// } catch (IOException e) {
			// e.printStackTrace();
			// }
		} else {
			if (!Parameter.grm_bin_flag) {
				heReader.XtX = new double[2][2];
				heReader.XtY = new double[2];
				// *************************************read grm file
				FileInputStream fin = null;
				try {
					fin = new FileInputStream(heReader.grmFile);
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

				BufferedReader grmFile = new BufferedReader(xover);

				try {
					HashMap<Double, Integer> cat = new HashMap<Double, Integer>();
					while ((line = grmFile.readLine()) != null) {
						String[] s = line.split(delim);
						int id1 = Integer.parseInt(s[0]) - 1;
						int id2 = Integer.parseInt(s[1]) - 1;
						if (id1 == id2)
							continue;
						if (!(heReader.flag[id1] & heReader.flag[id2]))
							continue;
						double ds = 0;
						if (heReader.heType[Parameter.he_sd]) {
							ds = (heReader.y[id1][1] - heReader.y[id2][1])
									* (heReader.y[id1][1] - heReader.y[id2][1]);
						} else if (heReader.heType[Parameter.he_ss]) {
							ds = (heReader.y[id1][1] + heReader.y[id2][1])
									* (heReader.y[id1][1] + heReader.y[id2][1]);
						} else if (heReader.heType[Parameter.he_cp]) {
							ds = heReader.y[id1][1] * heReader.y[id2][1];
						}

						if (cat.containsKey(heReader.y[id1][1])) {
							Integer I = (Integer) cat.get(heReader.y[id1][1]);
							I++;
							cat.put(heReader.y[id1][1], I);
						} else {
							cat.put(heReader.y[id1][1], 1);
						}
						if (cat.containsKey(heReader.y[id2][1])) {
							Integer I = (Integer) cat.get(heReader.y[id2][1]);
							I++;
							cat.put(heReader.y[id2][1], I);
						} else {
							cat.put(heReader.y[id2][1], 1);
						}
						heReader.yyProd += ds * ds;
						heReader.XtX[0][0]++;
						double g = Double.parseDouble(s[3]);
						heReader.XtX[0][1] += g;
						heReader.XtX[1][0] += g;
						heReader.XtX[1][1] += g * g;

						heReader.XtY[0] += ds;
						heReader.XtY[1] += g * ds;

						heReader.lambda.XYProd += g * ds;
						heReader.lambda.XXProd += g * g;
						heReader.lambda.SY += ds;
						heReader.lambda.SX += g;

						heReader.lambda.N++;
					}

					if (cat.size() == 2) {
						heReader.isCC = true;
						Set<Double> set = cat.keySet();
						Iterator<Double> it = set.iterator();
						Double k1 = it.next();
						Double k2 = it.next();
						Integer c1 = cat.get(k1);
						Integer c2 = cat.get(k2);

						if (k1 > k2) {
							heReader.P = c1.doubleValue()
									/ (c1.doubleValue() + c2.doubleValue());
						} else {
							heReader.P = c2.doubleValue()
									/ (c1.doubleValue() + c2.doubleValue());
						}
					}
					grmFile.close();

				} catch (IOException e) {
					e.printStackTrace();
				}
			} else {//grm.bin
				HashMap<Double, Integer> cat = new HashMap<Double, Integer>();

				FileInputStream fileStream = null;
				try {
					fileStream = new FileInputStream(Parameter.grm_bin);
				} catch (FileNotFoundException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				DataInputStream bigEndianDataStream = new DataInputStream(fileStream);
				LittleEndianDataInputStream littleEndianDataStream = new LittleEndianDataInputStream(bigEndianDataStream, Float.SIZE);

				for (int i = 0; i < heReader.flag.length; i++) {
					for (int j = 0; j <= i; j++) {
						float g = 0;
						try {
							if (littleEndianDataStream.available()>0) {
								g = littleEndianDataStream.readFloat();
							}
						} catch (IOException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}

						int id1 = i;
						int id2 = j;
						if (id1 == id2)
							continue;
						if (!(heReader.flag[id1] & heReader.flag[id2]))
							continue;
						double ds = 0;
						if (heReader.heType[Parameter.he_sd]) {
							ds = (heReader.y[id1][1] - heReader.y[id2][1])
									* (heReader.y[id1][1] - heReader.y[id2][1]);
						} else if (heReader.heType[Parameter.he_ss]) {
							ds = (heReader.y[id1][1] + heReader.y[id2][1])
									* (heReader.y[id1][1] + heReader.y[id2][1]);
						} else if (heReader.heType[Parameter.he_cp]) {
							ds = heReader.y[id1][1] * heReader.y[id2][1];
						}

						if (cat.containsKey(heReader.y[id1][1])) {
							Integer I = (Integer) cat.get(heReader.y[id1][1]);
							I++;
							cat.put(heReader.y[id1][1], I);
						} else {
							cat.put(heReader.y[id1][1], 1);
						}
						if (cat.containsKey(heReader.y[id2][1])) {
							Integer I = (Integer) cat.get(heReader.y[id2][1]);
							I++;
							cat.put(heReader.y[id2][1], I);
						} else {
							cat.put(heReader.y[id2][1], 1);
						}
						heReader.yyProd += ds * ds;
						heReader.XtX[0][0]++;
						heReader.XtX[0][1] += g;
						heReader.XtX[1][0] += g;
						heReader.XtX[1][1] += g * g;

						heReader.XtY[0] += ds;
						heReader.XtY[1] += g * ds;

						heReader.lambda.XYProd += g * ds;
						heReader.lambda.XXProd += g * g;
						heReader.lambda.SY += ds;
						heReader.lambda.SX += g;

						heReader.lambda.N++;

					}
				}
			}
		}

		for (int i = 0; i < heReader.XtX[0].length; i++) {
			for (int j = 0; j < i; j++) {
				heReader.XtX[i][j] = heReader.XtX[j][i];
			}
		}
	}

	public void Regression() {

		DecimalFormat fmt = new DecimalFormat("#.###E0");

		RealMatrix Mat_XtX = new Array2DRowRealMatrix(heReader.XtX);
		RealMatrix Mat_XtY = new Array2DRowRealMatrix(heReader.XtY);

		RealMatrix Mat_XtX_Inv = (new LUDecompositionImpl(Mat_XtX)).getSolver()
				.getInverse();
		Mat_B = Mat_XtX_Inv.multiply(Mat_XtY);
		RealMatrix Mat_K = Mat_B.transpose().multiply(Mat_XtY);

		double mse = (heReader.yyProd - Mat_K.getEntry(0, 0))
				/ (heReader.dim - heReader.mpheno.length - 1);
		RealMatrix v = Mat_XtX_Inv.scalarMultiply(mse);

		heReader.sb.append("HE mode: ");
		if (heReader.heType[Parameter.he_sd]) {
			heReader.sb.append("squared difference (yi-yj)^2\n");
		} else if (heReader.heType[Parameter.he_ss]) {
			heReader.sb.append("squared sum (yi+yj)^2\n");
		} else if (heReader.heType[Parameter.he_cp]) {
			heReader.sb.append("cross-product [yi-E(y)][yj-E(y)]\n");
		}
		heReader.sb.append("grm: " + heReader.grmFile + "\n");
		heReader.sb.append("grm id: " + heReader.grmID + "\n");
		heReader.sb.append("keep list: " + heReader.keepFile + "\n");
		heReader.sb.append("pheno: " + heReader.phenoFile + "\n");
		heReader.sb.append("mpheno: ");

		for (int i = 0; i < heReader.mpheno.length; i++) {
			heReader.sb.append(heReader.mpheno[i] + " ");
		}
		heReader.sb.append("\n");
		if (Parameter.INSTANCE.eh2Flag) {
			heReader.sb.append("Empirical h2: " + Parameter.INSTANCE.eh2 + "\n");
		}

		heReader.sb.append("reverse: " + heReader.reverse + "\n");
		if (heReader.k_button) {
			heReader.sb.append("k: " + heReader.k + "\n");
		}
		heReader.sb.append("\n========================\n");
		heReader.sb.append("Coef\t" + "Estimate \t" + "se" + "\n");

		for (int i = 0; i < Mat_B.getRowDimension(); i++) {
			heReader.sb.append("b" + i + "\t"
					+ fmt.format(Mat_B.getEntry(i, 0)) + "\t"
					+ fmt.format(Math.sqrt(v.getEntry(i, i))) + "\n");
		}

		if (!heReader.reverse) {
			double h_o = 0;
			if (heReader.heType[Parameter.he_sd]) {
				h_o = Mat_B.getEntry(1, 0) / Mat_B.getEntry(0, 0) * (-1);
			} else if (heReader.heType[Parameter.he_ss]) {
				h_o = Mat_B.getEntry(1, 0) / Mat_B.getEntry(0, 0);
			} else if (heReader.heType[Parameter.he_cp]) {
				h_o = Mat_B.getEntry(1, 0);
			}

			double u_b0 = Mat_B.getEntry(0, 0);
			double u_b1 = Mat_B.getEntry(1, 0);

			double v_b0 = v.getEntry(0, 0);
			double v_b1 = v.getEntry(1, 1);

			double v_ho = (u_b1 / u_b0)
					* (u_b1 / u_b0)
					* (v_b0 / (u_b0 * u_b0) + v_b1 / (u_b1 * u_b1) - 2
							* v.getEntry(0, 1) / (u_b0 * u_b1));
			if (heReader.heType[Parameter.he_cp]) {
				v_ho = Math.sqrt(v_b1);
			}
			heReader.sb.append("h2(o): " + fmt.format(h_o) + "\t"
					+ fmt.format(v_ho) + "\n");

			if (heReader.k_button && heReader.isCC) {
				NormalDistributionImpl Norm = new NormalDistributionImpl();
				double q = 0;
				try {
					q = Norm.inverseCumulativeProbability(1 - heReader.k);
				} catch (MathException e) {
					e.printStackTrace();
				}
				double z = 1 / (Math.sqrt(2 * 3.1416926))
						* Math.exp(-q * q / 2);
				double h_l = h_o * heReader.k * (1 - heReader.k) * heReader.k
						* (1 - heReader.k)
						/ (z * z * heReader.P * (1 - heReader.P));

				double f = (heReader.k * (1 - heReader.k) * heReader.k * (1 - heReader.k))
						/ (z * z * heReader.P * (1 - heReader.P));
				double v_hl = v_ho * f * f;
				heReader.sb.append("h2(l): " + fmt.format(h_l) + "\t"
						+ fmt.format(v_hl) + "\n");
			}
		}

		heReader.sb.append("\n");
		heReader.sb.append("variance-covariance\n");
		for (int i = 0; i < v.getRowDimension(); i++) {
			for (int j = 0; j <= i; j++) {
				heReader.sb.append(fmt.format(v.getEntry(i, j)) + " ");
			}
			heReader.sb.append("\n");
		}

		heReader.lambda.calLambda();
		heReader.sb
				.append("Decomposite b1 into the covariance and the variance:\n");
		heReader.sb.append("\ncov(Y, X): "
				+ fmt.format(heReader.lambda.getCov()) + "\n");
		heReader.sb.append("var(X): " + fmt.format(heReader.lambda.getVar())
				+ "\n");

		if (Parameter.INSTANCE.eh2Flag) {
			double Lmd = heReader.lambda.getLambda(-1 * Mat_B.getEntry(0, 0) * Parameter.INSTANCE.eh2);
			heReader.sb.append("Lambda: " + fmt.format(Lmd));
		}

		System.out.println(heReader.sb);

		if (heReader.output != null) {
			StringBuilder fsb = new StringBuilder();
			fsb.append(heReader.output);
			fsb.append(".he");
			File of = new File(fsb.toString());
			PrintWriter pw = null;
			try {
				pw = new PrintWriter(of);
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
			pw.append(heReader.sb);
			pw.close();
		}
	}

	public RealMatrix getB() {
		return Mat_B;
	}
}
