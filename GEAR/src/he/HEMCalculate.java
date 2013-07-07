package he;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Iterator;
import java.util.Set;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistributionImpl;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.Array2DRowRealMatrix;

import gear.CmdArgs;
import gear.ConstValues;
import gear.util.BinaryInputFile;
import gear.util.Logger;

public class HEMCalculate
{

	private final String delim = "\\s+";
	private HEMRead heMReader;
	private RealMatrix Mat_B;
	private int Len = 0;
	private double grmCutoff = 0;

	public HEMCalculate(HEMRead h)
	{
		heMReader = h;
		Calculate();
		
		for (BinaryInputFile grmBin : heMReader.binList)
		{
			grmBin.close();
		}
	}

	public void Calculate()
	{

		heMReader.lambda = new Lambda();

		if (CmdArgs.INSTANCE.getHEArgs().isAbsGrmCutoff())
		{
			grmCutoff = CmdArgs.INSTANCE.getHEArgs().AbsGrmCutoff();
		}
		else if (CmdArgs.INSTANCE.getHEArgs().isGrmCutoff())
		{
			grmCutoff = CmdArgs.INSTANCE.getHEArgs().GrmCutoff();
		}
		String line;

		Len = 0;
		for (int i = 0; i < heMReader.flag.length; i++)
			if (heMReader.flag[i])
				Len++;

		if (!heMReader.isBinGrm)
		{
			heMReader.XtX = new double[heMReader.grmList.size() + 1][heMReader.grmList
					.size() + 1];
			heMReader.XtY = new double[heMReader.grmList.size() + 1];

			try
			{
				int id1 = 0;
				int id2 = 0;
				for (int n = 0; n < heMReader.nRec; n++)
				{
					double[] row = new double[heMReader.grmList.size() + 1];
					row[0] = 1;
					boolean pass = false;
					for (int i = 0; i < heMReader.grmList.size(); i++)
					{
						BufferedReader grmFile_ = heMReader.grmList.get(i);
						if ((line = grmFile_.readLine()) == null)
						{
							continue;
						}

						String[] s = line.split(delim);
						if (i == 0)
						{
							id1 = Integer.parseInt(s[0]) - 1;
							id2 = Integer.parseInt(s[1]) - 1;
							if (id1 == id2)
							{
								pass = true;
							}
							if (!(heMReader.flag[id1] & heMReader.flag[id2]))
							{
								pass = true;
							}
						}
						row[1 + i] = Double.parseDouble(s[3]);
						if (heMReader.isSingleGrm)
						{
							if (CmdArgs.INSTANCE.getHEArgs().isAbsGrmCutoff())
							{
								if (Math.abs(row[1+i]) > grmCutoff)
								{
									pass = true;
								}
							}
							else if (CmdArgs.INSTANCE.getHEArgs().isGrmCutoff())
							{
								if (row[1+i] > grmCutoff)
								{
									pass = true;
								}
							}
						}
					}

					if (pass)
						continue;

					double ds = 0;

					switch (heMReader.heType)
					{
					case SD:
						ds = (heMReader.y[id1][1] - heMReader.y[id2][1]) * (heMReader.y[id1][1] - heMReader.y[id2][1]);
						break;
					case SS:
						ds = (heMReader.y[id1][1] + heMReader.y[id2][1]) * (heMReader.y[id1][1] + heMReader.y[id2][1]);
						break;
					case CP:
						ds = heMReader.y[id1][1] * heMReader.y[id2][1];
						break;
					default:
						break;
					}

					for (int j = 0; j < row.length; j++)
					{
						for (int k = 0; k <= j; k++)
						{
							heMReader.XtX[j][k] += row[j] * row[k];
						}
						heMReader.XtY[j] += ds * row[j];
					}
					heMReader.yyProd += ds * ds;
					heMReader.lambda.N++;
				}
				for (int i = 0; i < heMReader.grmList.size(); i++)
				{
					BufferedReader grmFile_ = heMReader.grmList.get(i);
					grmFile_.close();
				}
			}
			catch (IOException e)
			{
				if (CmdArgs.INSTANCE.getHEArgs().isGrmTxtList())
				{
					Logger.handleException(	e,
											"Error in reading Grm text List in Multiple regression for the HE model.");
				}
				if (CmdArgs.INSTANCE.getHEArgs().isGrmList())
				{
					Logger.handleException(	e,
											"Error in reading Grm gz List in Multiple regression for the HE model.");
				}
			}
		}
		else
		{
			heMReader.XtX = new double[heMReader.binList.size() + 1][heMReader.binList
					.size() + 1];
			heMReader.XtY = new double[heMReader.binList.size() + 1];

			double[] row = new double[heMReader.binList.size() + 1];
			row[0] = 1;
			for (int i = 0; i < heMReader.flag.length; i++)
			{
				for (int j = 0; j <= i; j++)
				{
					boolean pass = false;
					int id1 = i;
					int id2 = j;
					if (id1 == id2)
					{
						pass = true;
					}
					if (!(heMReader.flag[id1] & heMReader.flag[id2]))
					{
						pass = true;
					}
					for (int k = 0; k < heMReader.binList.size(); k++)
					{
						BinaryInputFile grmBin = heMReader.binList.get(0);

						if (grmBin.available() >= ConstValues.FLOAT_SIZE)
						{
							row[1 + k] = grmBin.readFloat();
						}

						if (heMReader.isSingleGrm)
						{
							if (CmdArgs.INSTANCE.getHEArgs().isAbsGrmCutoff())
							{
								if (Math.abs(row[1+i]) > grmCutoff)
								{
									pass = true;
								}
							}
							else if (CmdArgs.INSTANCE.getHEArgs().isGrmCutoff())
							{
								if (row[1+i] > grmCutoff)
								{
									pass = true;
								}
							}
						}
					}

					if (pass)
					{
						continue;
					}

					double ds = 0;

					switch (heMReader.heType)
					{
					case SD:
						ds = (heMReader.y[id1][1] - heMReader.y[id2][1]) * (heMReader.y[id1][1] - heMReader.y[id2][1]);
						break;
					case SS:
						ds = (heMReader.y[id1][1] + heMReader.y[id2][1]) * (heMReader.y[id1][1] + heMReader.y[id2][1]);
						break;
					case CP:
						ds = heMReader.y[id1][1] * heMReader.y[id2][1];
						break;
					default:
						// TODO: assert false or throw exception
						break;
					}

					for (int ii = 0; ii < row.length; ii++)
					{
						for (int jj = 0; jj <= ii; jj++)
						{
							heMReader.XtX[ii][jj] += row[ii] * row[jj];
						}
						heMReader.XtY[ii] += ds * row[ii];
					}
					heMReader.yyProd += ds * ds;
					heMReader.lambda.N++;
				}
			}
		}

		if (heMReader.cat.size() == 2)
		{
			heMReader.isCC = true;
			Set<Double> set = heMReader.cat.keySet();
			Iterator<Double> it = set.iterator();
			Double k1 = it.next();
			Double k2 = it.next();
			Integer c1 = heMReader.cat.get(k1);
			Integer c2 = heMReader.cat.get(k2);

			if (k1 > k2)
			{
				heMReader.P = c1.doubleValue() / (c1.doubleValue() + c2
						.doubleValue());
			}
			else
			{
				heMReader.P = c2.doubleValue() / (c1.doubleValue() + c2
						.doubleValue());
			}
		}

		for (int i = 0; i < heMReader.XtX.length; i++)
		{
			for (int j = 0; j < i; j++)
			{
				heMReader.XtX[j][i] = heMReader.XtX[i][j];
			}
		}

	}

	public void Regression()
	{

		DecimalFormat fmt = new DecimalFormat("#.###E0");

		RealMatrix Mat_XtX = new Array2DRowRealMatrix(heMReader.XtX);
		RealMatrix Mat_XtY = new Array2DRowRealMatrix(heMReader.XtY);

		RealMatrix Mat_XtX_Inv = (new LUDecompositionImpl(Mat_XtX)).getSolver()
				.getInverse();
		Mat_B = Mat_XtX_Inv.multiply(Mat_XtY);
		RealMatrix Mat_K = Mat_B.transpose().multiply(Mat_XtY);

		double mse = (heMReader.yyProd - Mat_K.getEntry(0, 0)) / (heMReader.lambda.N - heMReader.XtX.length);
		System.out.println("mse: " + mse);
		RealMatrix v = Mat_XtX_Inv.scalarMultiply(mse);

		heMReader.sb.append("HE mode: ");

		switch (heMReader.heType)
		{
		case SD:
			heMReader.sb.append("squared difference (yi-yj)^2\n");
			break;
		case SS:
			heMReader.sb.append("squared sum (yi+yj)^2\n");
			break;
		case CP:
			heMReader.sb.append("cross-product [yi-E(y)][yj-E(y)]\n");
			break;
		default:
			// TODO: assert false or throw exception
		}

		heMReader.sb.append("grm list file: " + heMReader.grmListFile + "\n");

		if (CmdArgs.INSTANCE.getHEArgs().isGrmCutoff())
		{
			heMReader.sb.append("grm cutoff is: " + grmCutoff + "\n");
		}

		if (CmdArgs.INSTANCE.getHEArgs().isAbsGrmCutoff())
		{
			heMReader.sb
					.append("grm absolute cutoff is: |" + grmCutoff + "|\n");
		}

		heMReader.sb.append("keep list: " + heMReader.keepFile + "\n");
		heMReader.sb.append("pheno: " + heMReader.phenoFile + "\n");
		heMReader.sb.append("mpheno: ");

		for (int i = 0; i < heMReader.mpheno.length; i++)
		{
			heMReader.sb.append(heMReader.mpheno[i] + " ");
		}
		heMReader.sb.append("\n");

		if (CmdArgs.INSTANCE.qcovar_file != null)
		{
			heMReader.sb
					.append("quantitative covariate file: " + CmdArgs.INSTANCE.qcovar_file + "\n");
			heMReader.sb.append("quantitative covariate index: ");
			if (CmdArgs.INSTANCE.qcovar_num == null)
			{
				heMReader.sb.append("all");
			}
			else
			{
				for (int i = 0; i < CmdArgs.INSTANCE.qcovar_num.length; i++)
				{
					heMReader.sb.append(CmdArgs.INSTANCE.qcovar_num[i] + " ");
				}
			}
			heMReader.sb.append("\n");
		}

		if (CmdArgs.INSTANCE.covar_file != null)
		{
			heMReader.sb
					.append("quality covariate file: " + CmdArgs.INSTANCE.covar_file + "\n");
			heMReader.sb.append("quality covariate index: ");
			if (CmdArgs.INSTANCE.covar_num == null)
			{
				heMReader.sb.append("all");
			}
			else
			{
				for (int i = 0; i < CmdArgs.INSTANCE.covar_num.length; i++)
				{
					heMReader.sb.append(CmdArgs.INSTANCE.covar_num[i] + " ");
				}
			}
			heMReader.sb.append("\n");
		}

		if (CmdArgs.INSTANCE.eh2Flag)
		{
			heMReader.sb.append("Empirical h2: " + CmdArgs.INSTANCE.eh2 + "\n");
		}

		heMReader.sb.append("reverse: " + heMReader.reverse + "\n");
		heMReader.sb.append("Scale : " + CmdArgs.INSTANCE.scale + "\n");
		if (heMReader.k_button)
		{
			heMReader.sb.append("k: " + heMReader.k + "\n");
		}

		heMReader.sb.append("In total " + Len + " matched individuals.\n");
		heMReader.sb
				.append("In total read " + ((int) (heMReader.lambda.N)) + " lines in the grm file.\n");
		heMReader.sb.append("\n========================\n");
		heMReader.sb.append("Coef\t" + "Estimate \t" + "se" + "\n");

		for (int i = 0; i < Mat_B.getRowDimension(); i++)
		{
			heMReader.sb.append("b" + i + "\t" + fmt.format(Mat_B
					.getEntry(i, 0)) + "\t" + fmt.format(Math.sqrt(v
					.getEntry(i, i))) + "\n");
		}

		if (!heMReader.reverse)
		{
			double[] h_o = new double[Mat_B.getRowDimension() - 1];
			double[] h_l = new double[Mat_B.getRowDimension() - 1];
			for (int i = 1; i < Mat_B.getRowDimension(); i++)
			{
				switch (heMReader.heType)
				{
				case SD:
					h_o[i - 1] = Mat_B.getEntry(i, 0) / Mat_B.getEntry(0, 0) * (-1);
					break;
				case SS:
					h_o[i - 1] = Mat_B.getEntry(i, 0) / Mat_B.getEntry(0, 0);
					break;
				case CP:
					h_o[i - 1] = Mat_B.getEntry(i, 0);
					break;
				default:
				}

				double u_b0 = Mat_B.getEntry(0, 0);
				double u_b1 = Mat_B.getEntry(i, 0);

				double v_b0 = v.getEntry(0, 0);
				double v_b1 = v.getEntry(i, i);

				double v_ho = (u_b1 / u_b0) * (u_b1 / u_b0) * (v_b0 / (u_b0 * u_b0) + v_b1 / (u_b1 * u_b1) - 2 * v
						.getEntry(0, i) / (u_b0 * u_b1));
				if (heMReader.heType == gear.HEType.CP)
				{
					v_ho = Math.sqrt(v_b1);
				}
				heMReader.sb.append("h2(o) for b" + i + ": " + fmt
						.format(h_o[i - 1]) + "\t" + fmt.format(v_ho) + "\n");

				if (heMReader.k_button && heMReader.isCC)
				{
					NormalDistributionImpl Norm = new NormalDistributionImpl();
					double q = 0;
					try
					{
						q = Norm.inverseCumulativeProbability(1 - heMReader.k);
					}
					catch (MathException e)
					{
						e.printStackTrace();
					}
					double z = 1 / (Math.sqrt(2 * 3.1416926)) * Math
							.exp(-q * q / 2);
					h_l[i - 1] = h_o[i - 1] * heMReader.k * (1 - heMReader.k) * heMReader.k * (1 - heMReader.k) / (z * z * heMReader.P * (1 - heMReader.P));

					double f = (heMReader.k * (1 - heMReader.k) * heMReader.k * (1 - heMReader.k)) / (z * z * heMReader.P * (1 - heMReader.P));
					double v_hl = v_ho * f * f;
					heMReader.sb
							.append("h2(l) for b" + i + ": " + fmt
									.format(h_l[i - 1]) + "\t" + fmt
									.format(v_hl) + "\n");
				}
			}
		}

		heMReader.sb.append("\n");
		heMReader.sb.append("variance-covariance\n");
		for (int i = 0; i < v.getRowDimension(); i++)
		{
			for (int j = 0; j <= i; j++)
			{
				heMReader.sb.append(fmt.format(v.getEntry(i, j)) + " ");
			}
			heMReader.sb.append("\n");
		}

		heMReader.lambda.calLambda();
		heMReader.sb
				.append("Decomposite b1 into the covariance and the variance:\n");
		heMReader.sb.append("\ncov(Y, X): " + fmt.format(heMReader.lambda
				.getCov()) + "\n");
		heMReader.sb
				.append("var(X): " + fmt.format(heMReader.lambda.getVar()) + "\n");

		if (CmdArgs.INSTANCE.eh2Flag)
		{
			double Lmd = heMReader.lambda
					.getLambda(-1 * Mat_B.getEntry(0, 0) * CmdArgs.INSTANCE.eh2);
			heMReader.sb.append("Lambda: " + fmt.format(Lmd));
		}

		Logger.printUserLog(heMReader.sb.toString());

		if (heMReader.output != null)
		{
			StringBuilder fsb = new StringBuilder();
			fsb.append(heMReader.output);
			fsb.append(".he");
			File of = new File(fsb.toString());
			PrintWriter pw = null;
			try
			{
				pw = new PrintWriter(of);
			}
			catch (FileNotFoundException e)
			{
				e.printStackTrace();
			}
			pw.append(heMReader.sb);
			pw.close();
		}
	}

	public RealMatrix getB()
	{
		return Mat_B;
	}
}
