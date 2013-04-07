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
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;
import java.util.logging.Level;
import java.util.zip.GZIPInputStream;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistributionImpl;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.Array2DRowRealMatrix;

import gear.Parameter;
import gear.util.FileProcessor;
import gear.util.Logger;
import he.endian.LittleEndianDataInputStream;


public class HECalculate {

	private final String delim = "\\s+";
	private HERead heReader;
	private RealMatrix Mat_B;
	private int Len = 0;
	private double grmCutoff = 0;

	public HECalculate(HERead h) {
		heReader = h;
		Calculate();
	}

	public void Calculate() {

		heReader.lambda = new Lambda();

		if (Parameter.INSTANCE.getHEParameter().isAbsGrmCutoff()) {
			grmCutoff = Parameter.INSTANCE.getHEParameter().AbsGrmCutoff();
		} else if (Parameter.INSTANCE.getHEParameter().isGrmCutoff()) {
			grmCutoff = Parameter.INSTANCE.getHEParameter().GrmCutoff();
		}
		String line;
		
		// ************************keep
		if (heReader.keepFile != null) {
			BufferedReader reader = FileProcessor.FileOpen(heReader.keepFile);
			boolean[] ff = new boolean[heReader.flag.length];
			Arrays.fill(ff, false);
			try {
				while ((line = reader.readLine()) != null) {
					String[] s = line.split(delim);
					StringBuilder sb = new StringBuilder(s[0] + "." + s[1]);
					if (heReader.ID2Idx.containsKey(sb.toString())) {
						int ii = heReader.ID2Idx.get(sb.toString());
						ff[ii] = true;
					}
				}
				reader.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
			for (int i = 0; i < ff.length; i++) {
				heReader.flag[i] &= ff[i];
			}
		}
		Len = 0;
		for (int i = 0; i < heReader.flag.length; i++)
			if (heReader.flag[i])
				Len++;
		
		// ************************************standardising
		double[] ss = new double[heReader.y[0].length - 1];
		double[] ssx = new double[heReader.y[0].length - 1];
		for (int i = 0; i < heReader.flag.length; i++) {
			if (!heReader.flag[i])
				continue;
			for (int j = 0; j < ss.length; j++) {
				ss[j] += heReader.y[i][j + 1];
				ssx[j] += heReader.y[i][j + 1] * heReader.y[i][j + 1];
			}
		}
		double[] sd = new double[ssx.length];
		for (int i = 0; i < sd.length; i++) {
			ss[i] /= Len;
			sd[i] = Math.sqrt((ssx[i] - Len * ss[i] * ss[i]) / (Len - 1));
		}
		if (Parameter.INSTANCE.scale) {
			Logger.printUserLog("Standardising phentoype.");
			for (int i = 0; i < heReader.flag.length; i++) {
				if (!heReader.flag[i])
					continue;
				for (int j = 1; j < heReader.y[i].length; j++) {
					heReader.y[i][j] = (heReader.y[i][j] - ss[j - 1]) / sd[j - 1];
				}
			}
		}
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
			if (Parameter.INSTANCE.getHEParameter().isGrmTxt()) {
				heReader.XtX = new double[2][2];
				heReader.XtY = new double[2];

				// *************************************read grm text file
				BufferedReader grmFile = FileProcessor.FileOpen(heReader.grmFile);

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
						double g = Double.parseDouble(s[3]);

						if (Parameter.INSTANCE.getHEParameter().isAbsGrmCutoff()) {
							if (Math.abs(g) > grmCutoff) continue;
						} else if (Parameter.INSTANCE.getHEParameter().isGrmCutoff()) {
							if (g > grmCutoff) continue;
						}
						double ds = 0;
						
						switch (heReader.heType) {
						case SD:
							ds = (heReader.y[id1][1] - heReader.y[id2][1]) *
							     (heReader.y[id1][1] - heReader.y[id2][1]);
							break;
						case SS:
							ds = (heReader.y[id1][1] + heReader.y[id2][1]) *
							     (heReader.y[id1][1] + heReader.y[id2][1]);
							break;
						case CP:
							ds = heReader.y[id1][1] * heReader.y[id2][1];
							break;
						default:
							// TODO: assert false or throw exception
							break;
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
			} else if (!Parameter.INSTANCE.getHEParameter().isGrmBinary()) {
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
						double g = Double.parseDouble(s[3]);

						if (Parameter.INSTANCE.getHEParameter().isAbsGrmCutoff()) {
							if (Math.abs(g) > grmCutoff) continue;
						} else if (Parameter.INSTANCE.getHEParameter().isGrmCutoff()) {
							if (g > grmCutoff) continue;
						}

						double ds = 0;
						
						switch (heReader.heType) {
						case SD:
							ds = (heReader.y[id1][1] - heReader.y[id2][1]) *
							     (heReader.y[id1][1] - heReader.y[id2][1]);
							break;
						case SS:
							ds = (heReader.y[id1][1] + heReader.y[id2][1]) *
							     (heReader.y[id1][1] + heReader.y[id2][1]);
							break;
						case CP:
							ds = heReader.y[id1][1] * heReader.y[id2][1];
							break;
						default:
							// TODO: assert false or throw exception
							break;
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
					fileStream = new FileInputStream(Parameter.INSTANCE.getHEParameter().getGrm());
				} catch (FileNotFoundException e) {
					Logger.printUserError("Cannot open the GRM file '" + Parameter.INSTANCE.getHEParameter().getGrm() + "'.");
					Logger.printUserError("Exception Message: " + e.getMessage());
					System.exit(1);
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
							Logger.printUserError("An exception occurred when reading the GRM file '" + Parameter.INSTANCE.getHEParameter().getGrm() + "'.");
							Logger.printUserError("Exception Message: " + e.getMessage());
							Logger.getDevLogger().log(Level.SEVERE, "Reading GRM file", e);
							System.exit(1);
						}

						int id1 = i;
						int id2 = j;
						if (id1 == id2)
							continue;
						if (!(heReader.flag[id1] & heReader.flag[id2]))
							continue;

						if (Parameter.INSTANCE.getHEParameter().isAbsGrmCutoff()) {
							if (Math.abs(g) > grmCutoff) continue;
						} else if (Parameter.INSTANCE.getHEParameter().isGrmCutoff()) {
							if (g > grmCutoff) continue;
						}

						double ds = 0;
						
						switch (heReader.heType) {
						case SD:
							ds = (heReader.y[id1][1] - heReader.y[id2][1]) *
							     (heReader.y[id1][1] - heReader.y[id2][1]);
							break;
						case SS:
							ds = (heReader.y[id1][1] + heReader.y[id2][1]) *
							     (heReader.y[id1][1] + heReader.y[id2][1]);
							break;
						case CP:
							ds = heReader.y[id1][1] * heReader.y[id2][1];
							break;
						default:
							// TODO: assert false or throw exception
							break;
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
		
		switch (heReader.heType) {
		case SD:
			heReader.sb.append("squared difference (yi-yj)^2\n");
			break;
		case SS:
			heReader.sb.append("squared sum (yi+yj)^2\n");
			break;
		case CP:
			heReader.sb.append("cross-product [yi-E(y)][yj-E(y)]\n");
			break;
		default:
			// TODO: assert false or throw exception
		}
		
		heReader.sb.append("grm: " + heReader.grmFile + "\n");
		heReader.sb.append("grm id: " + heReader.grmID + "\n");
		
		if (Parameter.INSTANCE.getHEParameter().isGrmCutoff()) {
			heReader.sb.append("grm cutoff is: " + grmCutoff + "\n");
		}

		if (Parameter.INSTANCE.getHEParameter().isAbsGrmCutoff()) {
			heReader.sb.append("grm absolute cutoff is: |" + grmCutoff + "|\n");
		}

		heReader.sb.append("keep list: " + heReader.keepFile + "\n");
		heReader.sb.append("pheno: " + heReader.phenoFile + "\n");
		heReader.sb.append("mpheno: ");

		for (int i = 0; i < heReader.mpheno.length; i++) {
			heReader.sb.append(heReader.mpheno[i] + " ");
		}
		heReader.sb.append("\n");

		if (Parameter.INSTANCE.qcovar_file != null) {
			heReader.sb.append("quantitative covariate file: " + Parameter.INSTANCE.qcovar_file + "\n");
			heReader.sb.append("quantitative covariate index: ");
			if (Parameter.INSTANCE.qcovar_num == null) {
				heReader.sb.append("all");
			} else {
				for (int i = 0; i < Parameter.INSTANCE.qcovar_num.length; i++) {
					heReader.sb.append(Parameter.INSTANCE.qcovar_num[i] + " ");
				}
			}
			heReader.sb.append("\n");
		}

		if (Parameter.INSTANCE.covar_file != null) {
			heReader.sb.append("quality covariate file: " + Parameter.INSTANCE.covar_file + "\n");
			heReader.sb.append("quality covariate index: ");
			if (Parameter.INSTANCE.covar_num == null) {
				heReader.sb.append("all");
			} else {
				for (int i = 0; i < Parameter.INSTANCE.covar_num.length; i++) {
					heReader.sb.append(Parameter.INSTANCE.covar_num[i] + " ");
				}
			}
			heReader.sb.append("\n");
		}

		if (Parameter.INSTANCE.eh2Flag) {
			heReader.sb.append("Empirical h2: " + Parameter.INSTANCE.eh2 + "\n");
		}

		heReader.sb.append("reverse: " + heReader.reverse + "\n");
		heReader.sb.append("Scale : " + Parameter.INSTANCE.scale + "\n");
		if (heReader.k_button) {
			heReader.sb.append("k: " + heReader.k + "\n");
		}
		
		heReader.sb.append("In total " + Len + " matched individuals.\n");
		heReader.sb.append("In total read " + ((int) (heReader.lambda.N)) + " lines in the grm file.\n");
		heReader.sb.append("\n========================\n");
		heReader.sb.append("Coef\t" + "Estimate \t" + "se" + "\n");

		for (int i = 0; i < Mat_B.getRowDimension(); i++) {
			heReader.sb.append("b" + i + "\t"
					+ fmt.format(Mat_B.getEntry(i, 0)) + "\t"
					+ fmt.format(Math.sqrt(v.getEntry(i, i))) + "\n");
		}

		if (!heReader.reverse) {
			double h_o = 0;
			
			switch (heReader.heType) {
			case SD:
				h_o = Mat_B.getEntry(1, 0) / Mat_B.getEntry(0, 0) * (-1);
				break;
			case SS:
				h_o = Mat_B.getEntry(1, 0) / Mat_B.getEntry(0, 0);
				break;
			case CP:
				h_o = Mat_B.getEntry(1, 0);
				break;
			default:
				// TODO: assert false or throw exception
			}
			
			double u_b0 = Mat_B.getEntry(0, 0);
			double u_b1 = Mat_B.getEntry(1, 0);

			double v_b0 = v.getEntry(0, 0);
			double v_b1 = v.getEntry(1, 1);

			double v_ho = (u_b1 / u_b0)
					* (u_b1 / u_b0)
					* (v_b0 / (u_b0 * u_b0) + v_b1 / (u_b1 * u_b1) - 2
							* v.getEntry(0, 1) / (u_b0 * u_b1));
			if (heReader.heType == gear.HEType.CP) {
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

		Logger.printUserLog(heReader.sb.toString());

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
