package he;

import gear.CmdArgs;
import gear.util.Logger;
import gear.util.NewIt;
import gear.util.Sample;
import he.endian.LittleEndianDataInputStream;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Random;
import java.util.Set;
import java.util.zip.GZIPInputStream;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistributionImpl;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;

public class HEPermutation
{

	private final String delim = "\\s+";
	private HERead heReader;

	private int[] index;
	private Random rnd = new Random(2012);
	private RealMatrix Mat_B;

	private double[] B;

	public HEPermutation(HERead h)
	{
		heReader = h;

		ArrayList<Integer> Idx = NewIt.newArrayList();
		for (int i = 0; i < heReader.flag.length; i++)
		{
			if (heReader.flag[i])
			{
				Idx.add(i);
			}
		}
		index = org.apache.commons.lang3.ArrayUtils.toPrimitive(Idx
				.toArray(new Integer[0]));
	}

	private double[][] permuteY()
	{
		double[][] y = new double[heReader.y.length][heReader.y[0].length];

		Sample.setSeed(rnd.nextLong());
		int[] idx = Sample.sample(index);

		int j = 0;
		for (int i = 0; i < y.length; i++)
		{
			if (!heReader.flag[i])
				continue;
			System.arraycopy(heReader.y[i], 0, y[idx[j]], 0, y[idx[j]].length);
			j++;
		}

		return y;
	}

	public void Permutation()
	{
		DecimalFormat fmt = new DecimalFormat("#.###E0");

		B = new double[heReader.perm];
		if (heReader.permFlag)
		{
			for (int i = 0; i < heReader.perm; i++)
			{
				this.Calculate();
				Regression();
				B[i] = Mat_B.getEntry(1, 0);
				Logger.printUserLog("Permutation " + i + ": b1 = " + B[i]);
			}
		}

		Arrays.sort(B);
		for (int i = 0; i < B.length; i++)
		{
			Logger.printUserLog(Double.toString(B[i]));
		}

		DescriptiveStatistics ds = new DescriptiveStatistics(B);
		Logger.printUserLog("Mean : " + fmt.format(ds.getMean()));
		Logger.printUserLog("Standard deviation : "
				+ fmt.format(ds.getStandardDeviation()));

	}

	public void Calculate()
	{
		String line;

		heReader.lambda = new Lambda();
		if (heReader.reverse)
		{
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
		} else
		{
			if (!CmdArgs.INSTANCE.getHEArgs().isGrmBinary())
			{
				heReader.XtX = new double[2][2];
				heReader.XtY = new double[2];
				// *************************************read grm file
				FileInputStream fin = null;
				try
				{
					fin = new FileInputStream(heReader.grmFile);
				} catch (FileNotFoundException e1)
				{
					e1.printStackTrace();
				}
				GZIPInputStream gzis = null;
				try
				{
					gzis = new GZIPInputStream(fin);
				} catch (IOException e1)
				{
					e1.printStackTrace();
				}
				InputStreamReader xover = new InputStreamReader(gzis);

				BufferedReader grmFile = new BufferedReader(xover);

				try
				{
					HashMap<Double, Integer> cat = new HashMap<Double, Integer>();
					while ((line = grmFile.readLine()) != null)
					{
						String[] s = line.split(delim);
						int id1 = Integer.parseInt(s[0]) - 1;
						int id2 = Integer.parseInt(s[1]) - 1;
						if (id1 == id2)
							continue;
						if (!(heReader.flag[id1] & heReader.flag[id2]))
							continue;
						double ds = 0;

						switch (heReader.heType)
						{
						case SD:
							ds = (heReader.y[id1][1] - heReader.y[id2][1])
									* (heReader.y[id1][1] - heReader.y[id2][1]);
							break;
						case SS:
							ds = (heReader.y[id1][1] + heReader.y[id2][1])
									* (heReader.y[id1][1] + heReader.y[id2][1]);
							break;
						case CP:
							ds = heReader.y[id1][1] * heReader.y[id2][1];
							break;
						default:
							// TODO: assert false or throw exception
						}

						if (cat.containsKey(heReader.y[id1][1]))
						{
							Integer I = (Integer) cat.get(heReader.y[id1][1]);
							I++;
							cat.put(heReader.y[id1][1], I);
						} else
						{
							cat.put(heReader.y[id1][1], 1);
						}
						if (cat.containsKey(heReader.y[id2][1]))
						{
							Integer I = (Integer) cat.get(heReader.y[id2][1]);
							I++;
							cat.put(heReader.y[id2][1], I);
						} else
						{
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

					if (cat.size() == 2)
					{
						heReader.isCC = true;
						Set<Double> set = cat.keySet();
						Iterator<Double> it = set.iterator();
						Double k1 = it.next();
						Double k2 = it.next();
						Integer c1 = cat.get(k1);
						Integer c2 = cat.get(k2);

						if (k1 > k2)
						{
							heReader.P = c1.doubleValue()
									/ (c1.doubleValue() + c2.doubleValue());
						} else
						{
							heReader.P = c2.doubleValue()
									/ (c1.doubleValue() + c2.doubleValue());
						}
					}
					grmFile.close();

				} catch (IOException e)
				{
					e.printStackTrace();
				}
			} else
			{// grm.bin
				HashMap<Double, Integer> cat = new HashMap<Double, Integer>();

				FileInputStream fileStream = null;
				try
				{
					fileStream = new FileInputStream(CmdArgs.INSTANCE
							.getHEArgs().getGrm());
				} catch (FileNotFoundException e)
				{
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				DataInputStream bigEndianDataStream = new DataInputStream(
						fileStream);
				LittleEndianDataInputStream littleEndianDataStream = new LittleEndianDataInputStream(
						bigEndianDataStream, Float.SIZE);

				for (int i = 0; i < heReader.flag.length; i++)
				{
					for (int j = 0; j <= i; j++)
					{
						float g = 0;
						try
						{
							if (littleEndianDataStream.available() > 0)
							{
								g = littleEndianDataStream.readFloat();
							}
						} catch (IOException e)
						{
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

						switch (heReader.heType)
						{
						case SD:
							ds = (heReader.y[id1][1] - heReader.y[id2][1])
									* (heReader.y[id1][1] - heReader.y[id2][1]);
							break;
						case SS:
							ds = (heReader.y[id1][1] + heReader.y[id2][1])
									* (heReader.y[id1][1] + heReader.y[id2][1]);
							break;
						case CP:
							ds = heReader.y[id1][1] * heReader.y[id2][1];
							break;
						default:
							// TODO: assert false or throw exception
						}

						if (cat.containsKey(heReader.y[id1][1]))
						{
							Integer I = (Integer) cat.get(heReader.y[id1][1]);
							I++;
							cat.put(heReader.y[id1][1], I);
						} else
						{
							cat.put(heReader.y[id1][1], 1);
						}
						if (cat.containsKey(heReader.y[id2][1]))
						{
							Integer I = (Integer) cat.get(heReader.y[id2][1]);
							I++;
							cat.put(heReader.y[id2][1], I);
						} else
						{
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

		for (int i = 0; i < heReader.XtX[0].length; i++)
		{
			for (int j = 0; j < i; j++)
			{
				heReader.XtX[i][j] = heReader.XtX[j][i];
			}
		}
	}

	public void Regression()
	{
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

		for (int i = 0; i < Mat_B.getRowDimension(); i++)
		{
			heReader.sb.append("b" + i + "\t"
					+ fmt.format(Mat_B.getEntry(i, 0)) + "\t"
					+ fmt.format(Math.sqrt(v.getEntry(i, i))) + "\n");
		}

		if (!heReader.reverse)
		{
			double h_o = 0;

			switch (heReader.heType)
			{
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
			if (heReader.heType == gear.HEType.CP)
			{
				v_ho = Math.sqrt(v_b1);
			}
			heReader.sb.append("h2(o): " + fmt.format(h_o) + "\t"
					+ fmt.format(v_ho) + "\n");

			if (heReader.k_button && heReader.isCC)
			{
				NormalDistributionImpl Norm = new NormalDistributionImpl();
				double q = 0;
				try
				{
					q = Norm.inverseCumulativeProbability(1 - heReader.k);
				} catch (MathException e)
				{
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
		for (int i = 0; i < v.getRowDimension(); i++)
		{
			for (int j = 0; j <= i; j++)
			{
				heReader.sb.append(fmt.format(v.getEntry(i, j)) + " ");
			}
			heReader.sb.append("\n");
		}
		// System.out.println(heReader.sb);

		if (heReader.output != null)
		{
			StringBuilder fsb = new StringBuilder();
			fsb.append(heReader.output);
			fsb.append(".he");
			File of = new File(fsb.toString());
			PrintWriter pw = null;
			try
			{
				pw = new PrintWriter(of);
			} catch (FileNotFoundException e)
			{
				e.printStackTrace();
			}
			pw.append(heReader.sb);
			pw.close();
		}
	}

}
