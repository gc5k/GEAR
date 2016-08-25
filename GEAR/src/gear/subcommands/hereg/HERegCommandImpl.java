package gear.subcommands.hereg;

import java.io.PrintStream;
import java.text.DecimalFormat;

import org.apache.commons.math.stat.StatUtils;
import org.apache.commons.math.stat.regression.OLSMultipleLinearRegression;
import org.apache.commons.math.stat.regression.SimpleRegression;

import gear.ConstValues;
import gear.data.InputDataSet2;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.BinaryInputFile;
import gear.util.BufferedReader;
import gear.util.FileUtil;
import gear.util.Logger;

public class HERegCommandImpl extends CommandImpl {

	@Override
	public void execute(CommandArguments cmdArgs)
	{
		heArgs = (HERegCommandArguments) cmdArgs;

		data = new InputDataSet2();
		int fidx = 1;

		if (heArgs.isGRMList())
		{
			data.addFile(heArgs.getGrmList()[0]+".grm.id");
			data.addFile(heArgs.getPhenotypeFile(), heArgs.getPhenotypeIdx());

			if (heArgs.getCovFile() != null)
			{
				data.addFile(heArgs.getCovFile(), heArgs.getCovNumber());
				fidx++;
				covFileIdx = fidx;
			}
			if (heArgs.getKeepFile() != null)
			{
				data.addFile(heArgs.getKeepFile());
				fidx++;
				keepFileIdx = fidx;
			}
			data.LineUpFiles();

			Y = readPhenotype();
			A3 = readGrmList(data.getFileSampleSize(grmFileIdx));
			Logger.printUserLog(A3.length + " variance components included.");
			Logger.printUserLog("");

			if ( heArgs.getCovFile() != null)
			{
				X = readCovar();
				adjustY();
			}
			
			if (heArgs.isScale())
			{
				Y = StatUtils.normalize(Y);
			}
			HEMultple();
		}
		else
		{
			data.addFile(heArgs.getGrmID());
			data.addFile(heArgs.getPhenotypeFile(), heArgs.getPhenotypeIdx());

			if (heArgs.getCovFile() != null)
			{
				data.addFile(heArgs.getCovFile(), heArgs.getCovNumber());
				fidx++;
				covFileIdx = fidx;
			}
			if (heArgs.getKeepFile() != null)
			{
				data.addFile(heArgs.getKeepFile());
				fidx++;
				keepFileIdx = fidx;
			}
			data.LineUpFiles();

			Y = readPhenotype();

			readGrm(data.getFileSampleSize(grmFileIdx));
			Logger.printUserLog("1 variance component included.");
			Logger.printUserLog("");

			if (heArgs.getCovFile() != null)
			{
				X = readCovar();
				adjustY();
			}

			if (heArgs.isScale())
			{
				Y = StatUtils.normalize(Y);
			}
			HESingle();
		}
		
		SummaryA3();
		
		printout();
	}

	private void printout()
	{
		String Fout = heArgs.getOutRoot()  + ".he";
		PrintStream HE = FileUtil.CreatePrintStream(Fout);
		HE.println();
		HE.println("GRM\tEffective samples\tEffective markers");
		HE.println("-------------------------------------------------------------------");

		Logger.printUserLog("");
		Logger.printUserLog("GRM\tEffective samples\tEffective markers");
		Logger.printUserLog("-------------------------------------------------------------------");

		double mA_= 0;
		double vA_= 0;
		for (int i = 0; i < mA.length; i++)
		{
			HE.println("GRM" + (i+1) + "\t" + fmt.format(-1/mA[i]) + "\t" + fmt.format(1/vA[i][i]));
			Logger.printUserLog("GRM" + (i+1) + "\t" + fmt.format(-1/mA[i]) + "\t" + fmt.format(1/vA[i][i]));

			mA_ += -1/mA[i];
			vA_ = 1/vA[i][i];
		}
		
		if (mA.length > 1)
		{
			HE.println("Ave.\t" + fmt.format(mA_) + "\t" + fmt.format(vA_));
			Logger.printUserLog("Ave.\t" + fmt.format(mA_) + "\t" + fmt.format(vA_));
		}

		HE.println("");
		HE.println("Covariance structure for GRM");
		HE.println("-------------------------------------------------------------------");
		
		Logger.printUserLog("");
		Logger.printUserLog("Covariance structure for GRM");
		Logger.printUserLog("-------------------------------------------------------------------");

		for (int i = 0; i < vA.length; i++)
		{
			String s = new String();
			for (int j = 0; j <= i; j++)
			{
				s += fmtE.format(vA[i][j]) + "\t";
			}
			HE.println(s);
			Logger.printUserLog(s);
		}

		HE.println("-------------------------------------------------------------------");
		HE.println("");
		HE.println("-------------------------------------------------------------------");
		HE.println("HE regression");
		HE.println("Para\tEstimate\tSE\tPROP");
		HE.println("-------------------------------------------------------------------");

		Logger.printUserLog("-------------------------------------------------------------------");
		Logger.printUserLog("");
		Logger.printUserLog("-------------------------------------------------------------------");
		Logger.printUserLog("HE regression");
		Logger.printUserLog("Para\tEstimate\tSE\tPROP");
		Logger.printUserLog("-------------------------------------------------------------------");

		for (int i = 0; i < beta.length; i++)
		{
			if (i == 0)
			{
				HE.println("Mean\t" + fmt.format(beta[i]) + "\t" + fmt.format( Math.sqrt(beta_v[i][i])));
				Logger.printUserLog("Mean\t" + fmt.format(beta[i]) + "\t" + fmt.format( Math.sqrt(beta_v[i][i])));
			}
			else
			{
				HE.println("Beta" + i +"\t" + fmt.format(beta[i]) + "\t" + fmt.format(Math.sqrt(beta_v[i][i])) + "\t" + fmt.format(-1*beta[i]/beta[0]));
				Logger.printUserLog("Beta" + i +"\t" + fmt.format(beta[i]) + "\t" + fmt.format(Math.sqrt(beta_v[i][i])) + "\t" + fmt.format(-1*beta[i]/beta[0]));
			}
		}
		HE.println("-------------------------------------------------------------------");
		Logger.printUserLog("-------------------------------------------------------------------");

		HE.println("Covariance structure for HE");
		Logger.printUserLog("Covariance structure for HE");

		for (int i = 0; i < beta_v.length; i++)
		{
			String s = new String();
			for (int j = 0; j <= i; j++)
			{
				s += fmtE.format(beta_v[i][j]) + "\t";
			}
			HE.println(s);
			Logger.printUserLog(s);
		}
		HE.close();
	}
	
	private void HEMultple()
	{
		OLSMultipleLinearRegression mod = new OLSMultipleLinearRegression();

		double[] yy = new double[Y.length * (Y.length-1)/2];
		double[][] xx = new double[Y.length * (Y.length -1)/2][A3.length];
		int cnt = 0;
		for (int i = 0; i < Y.length; i++)
		{
			for (int j = 0; j < i; j++)
			{
				if (heArgs.isSD())
				{
					yy[cnt] = (Y[i] - Y[j]) * (Y[i] - Y[j]);
				}
				else if (heArgs.isCP())
				{
					yy[cnt] = Y[i]*Y[j];
				}
				else
				{
					yy[cnt] = (Y[i] + Y[j]) * (Y[i] + Y[j]);
				}
				for (int k = 0; k < A3.length; k++)
				{
					xx[cnt][k] = A3[k][i][j];				
				}
				cnt++;
			}
		}

		mod.newSampleData(yy, xx);
		beta = mod.estimateRegressionParameters();
		beta_v = mod.estimateRegressionParametersVariance();
	}

	private void SummaryA3()
	{
		mA = new double[A3.length];
		vA = new double[A3.length][A3.length];
		for (int i = 0; i < A3.length; i++)
		{
			for (int j = 0; j <=i; j++)
			{
				double[] a = calA2(A3[i], A3[j]);
				if (i == j) mA[j] = a[0];
				vA[i][j] = vA[j][i] = a[1];
			}
		}
	}

	private double[] calA2(double[][] b1, double[][] b2)
	{
		double Sx1 = 0;
		double Sx2 = 0;
		double Sxx = 0;

		for (int i = 1; i < b1.length; i++)
		{
			for (int j = 0; j < i; j++)
			{
				Sx1 += b1[i][j];
				Sx2 += b2[i][j];
				Sxx += b1[i][j] * b2[i][j];
			}
		}
		int N = b1.length * (b1.length -1)/2;

		double[] br = {Sx1/N, (Sxx/N - Sx1/N*Sx2/N)};
		return br;
	}

	private void HESingle()
	{
		SimpleRegression sReg = new SimpleRegression();

		int cnt = 0;
		for (int i = 0; i < Y.length; i++)
		{
			for (int j = 0; j < i; j++)
			{
				double yy = 0;
				if (heArgs.isSD())
				{
					yy = (Y[i] - Y[j]) * (Y[i] - Y[j]);
				}
				else if (heArgs.isCP())
				{
					yy = Y[i]*Y[j];
				}
				else
				{
					yy = (Y[i] + Y[j]) * (Y[i] + Y[j]);
				}
				if (heArgs.isGRMcut() && A3[0][i][j] > heArgs.getGRMcutoff()) continue;
				sReg.addData(A3[0][i][j], yy);
				cnt++;
			}
		}

		Logger.printUserLog(cnt + " observations for HE regression.");
		beta = new double[2];
		beta_v = new double[2][2];

		beta[0] = sReg.getIntercept();
		beta[1] = sReg.getSlope();
		beta_v[0][0] = sReg.getInterceptStdErr();
		beta_v[1][1] = sReg.getSlopeStdErr();
	}
	
	
	private double[][] readCovar()
	{
		int[] covIdx = data.getMatchedSubjectIdx(covFileIdx);
		double[][] x = new double[covIdx.length][heArgs.getCovNumber().length];
		for (int subjectIdx = 0; subjectIdx < covIdx.length; subjectIdx++)
		{
			x[subjectIdx] = data.getVariable(covFileIdx, covIdx[subjectIdx], heArgs.getCovNumber());
		}
		return x;
	}

	private double[] readPhenotype()
	{
		int[] pheIdx = data.getMatchedSubjectIdx(pheFileIdx);
		double[] y = new double[pheIdx.length];
		for(int subjectIdx = 0; subjectIdx < y.length; subjectIdx++)
		{
			y[subjectIdx] = data.getVariable(pheFileIdx, pheIdx[subjectIdx], heArgs.getPhenotypeIdx()[0]);
		}
		return y;
	}

	private void adjustY()
	{		
		OLSMultipleLinearRegression mod = new OLSMultipleLinearRegression();
		mod.newSampleData(Y, X);
		double[] res = mod.estimateResiduals();
		System.arraycopy(res, 0, Y, 0, Y.length);
	}

	private double[][][] lineUpGenotype(double[][][] B)
	{
		int[] subIdx = data.getMatchedSubjectIdx(grmFileIdx);
		double[][][] b = new double[B.length][subIdx.length][subIdx.length];
		
		for (int i = 0; i < b.length; i++)
		{
			b[i] = lineUpMatrix(B[i], subIdx);
		}
		return b;
	}

	private double[][] lineUpMatrix(double[][] B, int[] subIdx)
	{
		double[][] a = new double[subIdx.length][subIdx.length];

		for (int i = 0; i < subIdx.length; i++)
		{
			for (int j = 0; j < subIdx.length; j++)
			{
				a[i][j] = a[j][i] = B[subIdx[i]][subIdx[j]];
			}
		}
		return a;
	}

	private double[][][] readGrmList(int numSubjects)
	{
		double[][][] B = new double[heArgs.getGrmList().length][numSubjects][numSubjects];
		String[] tokens = null;
		for (int z = 0; z < heArgs.getGrmList().length; z++)
		{
			BufferedReader reader = BufferedReader.openGZipFile(heArgs.getGrmList()[z] + ".grm.gz", "GRM (gzip)");

			for (int i = 0; i < B[z].length; i++)
			{
				for (int j = 0; j <= i; j++)
				{
					if ((tokens = reader.readTokens(4)) != null) 
					{
						B[z][i][j] = B[z][j][i] = Double.parseDouble(tokens[3]);
					}
				}
			}
			reader.close();
		}
		return lineUpGenotype(B);
	}

	private void readGrm(int numSubjects)
	{
		if (heArgs.getGrmBin() != null)
		{
			readGrmBin(heArgs.getGrmBin(), numSubjects);
		}
		else
		{
			BufferedReader reader = heArgs.getGrmText() == null ?
					BufferedReader.openGZipFile(heArgs.getGrmGZ(), "GRM (gzip)") :
					BufferedReader.openTextFile(heArgs.getGrmText(), "GRM");
			readGrm(reader, numSubjects);
		}
	}

	private void readGrmBin(String fileName, int numSubjects)
	{
		BinaryInputFile grmBin = new BinaryInputFile(fileName, "GRM (binary)", /*littleEndian*/true);
		double[][] B = new double[numSubjects][numSubjects];
		Logger.printUserLog("Constructing A matrix: a " + numSubjects + " X " + numSubjects + " matrix.");
		for (int i = 0; i < B.length; i++) 
		{
			for (int j = 0; j <= i; j++)
			{
				if (grmBin.available() >= ConstValues.FLOAT_SIZE)
				{
					B[i][j] = B[j][i] = grmBin.readFloat();
				}
			}
		}
		grmBin.close();
		A3 = new double[1][data.getMatchedSubjectIdx(grmFileIdx).length][data.getMatchedSubjectIdx(grmFileIdx).length];
		A3[0] = lineUpMatrix(B, data.getMatchedSubjectIdx(grmFileIdx));
	}

	private void readGrm(BufferedReader reader, int numSubjects)
	{
		double[][] B = new double[numSubjects][numSubjects];
		String[] tokens = null;
		for (int i = 0; i < B.length; i++)
		{
			for (int j = 0; j <= i; j++)
			{
				if ((tokens = reader.readTokens(4)) != null) 
				{
					B[i][j] = B[j][i] = Double.parseDouble(tokens[3]);
				}
			}
		}
		reader.close();
		A3 = new double[1][data.getMatchedSubjectIdx(grmFileIdx).length][data.getMatchedSubjectIdx(grmFileIdx).length];
		A3[0] = lineUpMatrix(B, data.getMatchedSubjectIdx(grmFileIdx));
	}

	private int grmFileIdx = 0;
	private int pheFileIdx = 1;
	private int covFileIdx = 2;
	private int keepFileIdx = 3;

	private double[][][] A3;
	
	private double[] Y = null;
	private double[][] X = null;
	private double[] beta = null;
	private double[][] beta_v = null;

	private HERegCommandArguments heArgs = null;
	private InputDataSet2 data = null;
	
	private double[] mA;
	private double[][] vA;
	
	private DecimalFormat fmt = new DecimalFormat("0.0000");
	private DecimalFormat fmtE = new DecimalFormat("0.00E0");
	
}
