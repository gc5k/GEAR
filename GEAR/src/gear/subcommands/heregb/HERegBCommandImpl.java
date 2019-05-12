package gear.subcommands.heregb;

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

public class HERegBCommandImpl extends CommandImpl {

	@Override
	public void execute(CommandArguments cmdArgs) {
		heBArgs = (HERegBCommandArguments) cmdArgs;

		data = new InputDataSet2();

		if (heBArgs.isGRMList()) {
			data.addFile(heBArgs.getGrmList()[0] + ".grm.id");
			data.addFile(heBArgs.getPhenotypeFile(), heBArgs.getPhenotypeIdx());

			if (heBArgs.getCovFile() != null) {
				data.addFile(heBArgs.getCovFile(), heBArgs.getCovNumber());
			}
			if (heBArgs.getKeepFile() != null) {
				data.addFile(heBArgs.getKeepFile());
			}
			data.LineUpFiles();

			Y = readPhenotypes();
			A3 = readGrmList(data.getFileSampleSize(grmFileIdx));
			Logger.printUserLog(A3.length + " variance components included.");
			Logger.printUserLog("");

			if (heBArgs.getCovFile() != null) {
				X = readCovar();
				adjustY();
			}

			if (heBArgs.isScale()) {
				Y[0] = StatUtils.normalize(Y[0]);
				Y[1] = StatUtils.normalize(Y[1]);
			}
			HEMultple();
		} else {
			data.addFile(heBArgs.getGrmID());
			data.addFile(heBArgs.getPhenotypeFile(), heBArgs.getPhenotypeIdx());

			if (heBArgs.getCovFile() != null) {
				data.addFile(heBArgs.getCovFile(), heBArgs.getCovNumber());
			}
			if (heBArgs.getKeepFile() != null) {
				data.addFile(heBArgs.getKeepFile());
			}
			data.LineUpFiles();

			Y = readPhenotypes();

			readGrm(data.getFileSampleSize(grmFileIdx));
			Logger.printUserLog("1 variance component included.");
			Logger.printUserLog("");

			if (heBArgs.getCovFile() != null) {
				X = readCovar();
				adjustY();
			}

			if (heBArgs.isScale()) {
				double[] y0 = new double[Y[0].length], y1 = new double[Y[1].length];
				System.arraycopy(Y[0], 0, y0, 0, Y[0].length);
				System.arraycopy(Y[1], 0, y1, 0, Y[1].length);

				y0 = StatUtils.normalize(y0);
				y1 = StatUtils.normalize(y1);
				System.arraycopy(y0, 0, Y[0], 0, y0.length);
				System.arraycopy(y1, 0, Y[1], 0, y1.length);

			}
			HESingle();

			if (heBArgs.isJackknife()) {
				HESingleJackknife();
			}
		}

		SummaryA3();
		printout();
	}

	private void printout() {
		String Fout = heBArgs.getOutRoot() + ".he";
		PrintStream HE = FileUtil.CreatePrintStream(Fout);
		HE.println();
		HE.println("GRM\tEffective samples\tEffective markers");
		HE.println("-------------------------------------------------------------------");

		Logger.printUserLog("");
		Logger.printUserLog("GRM\tEffective samples\tEffective markers");
		Logger.printUserLog("-------------------------------------------------------------------");

		double mA_ = 0;
		double vA_ = 0;
		for (int i = 0; i < mA.length; i++) {
			HE.println("GRM" + (i + 1) + "\t" + fmt.format(-1 / mA[i]) + "\t" + fmt.format(1 / vA[i][i]));
			Logger.printUserLog("GRM" + (i + 1) + "\t" + fmt.format(-1 / mA[i]) + "\t" + fmt.format(1 / vA[i][i]));

			mA_ += -1 / mA[i];
			vA_ = 1 / vA[i][i];
		}

		if (mA.length > 1) {
			HE.println("Ave.\t" + fmt.format(mA_) + "\t" + fmt.format(vA_));
			Logger.printUserLog("Ave.\t" + fmt.format(mA_) + "\t" + fmt.format(vA_));
		}

		HE.println("");
		HE.println("Covariance structure for GRM");
		HE.println("-------------------------------------------------------------------");

		Logger.printUserLog("");
		Logger.printUserLog("Covariance structure for GRM");
		Logger.printUserLog("-------------------------------------------------------------------");

		for (int i = 0; i < vA.length; i++) {
			String s = new String();
			for (int j = 0; j <= i; j++) {
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

		for (int i = 0; i < beta.length; i++) {
			if (i == 0) {
				HE.println("Mean\t" + fmt.format(beta[i]) + "\t" + fmt.format(Math.sqrt(beta_v[i][i])));
				Logger.printUserLog("Mean\t" + fmt.format(beta[i]) + "\t" + fmt.format(Math.sqrt(beta_v[i][i])));
			} else {
				HE.println("Beta" + i + "\t" + fmt.format(beta[i]) + "\t" + fmt.format(Math.sqrt(beta_v[i][i])) + "\t"
						+ fmt.format(-1 * beta[i] / beta[0]));
				Logger.printUserLog("Beta" + i + "\t" + fmt.format(beta[i]) + "\t" + fmt.format(Math.sqrt(beta_v[i][i]))
						+ "\t" + fmt.format(-1 * beta[i] / beta[0]));
			}
		}
		HE.println("-------------------------------------------------------------------");
		Logger.printUserLog("-------------------------------------------------------------------");

		HE.println("Covariance structure for HE");
		Logger.printUserLog("Covariance structure for HE");

		for (int i = 0; i < beta_v.length; i++) {
			String s = new String();
			for (int j = 0; j <= i; j++) {
				s += fmtE.format(beta_v[i][j]) + "\t";
			}
			HE.println(s);
			Logger.printUserLog(s);
		}
		HE.close();
	}

	private void HEMultple() {
		OLSMultipleLinearRegression mod = new OLSMultipleLinearRegression();

		double[] yy = new double[Y[0].length * (Y[0].length - 1) / 2];
		double[][] xx = new double[Y[0].length * (Y[0].length - 1) / 2][A3.length];
		int cnt = 0;
		for (int i = 0; i < Y[0].length; i++) {
			for (int j = 0; j < i; j++) {
				if (heBArgs.isSD()) {
					yy[cnt] = (Y[0][i] - Y[0][j]) * (Y[1][i] - Y[1][j]);
				} else if (heBArgs.isCP()) {
					yy[cnt] = (Y[0][i] * Y[1][j] + Y[1][i] * Y[0][j]) / 2;
				} else {
					yy[cnt] = (Y[0][i] + Y[1][j]) * (Y[1][i] + Y[0][j]);
				}
				for (int k = 0; k < A3.length; k++) {
					xx[cnt][k] = A3[k][i][j];
				}
				cnt++;
			}
		}

		mod.newSampleData(yy, xx);
		beta = mod.estimateRegressionParameters();
		beta_v = mod.estimateRegressionParametersVariance();
	}

	private void SummaryA3() {
		mA = new double[A3.length];
		vA = new double[A3.length][A3.length];
		for (int i = 0; i < A3.length; i++) {
			for (int j = 0; j <= i; j++) {
				double[] a = calA2(A3[i], A3[j]);
				if (i == j)
					mA[j] = a[0];
				vA[i][j] = vA[j][i] = a[1];
			}
		}
	}

	private double[] calA2(double[][] b1, double[][] b2) {
		double Sx1 = 0;
		double Sx2 = 0;
		double Sxx = 0;

		for (int i = 1; i < b1.length; i++) {
			for (int j = 0; j < i; j++) {
				Sx1 += b1[i][j];
				Sx2 += b2[i][j];
				Sxx += b1[i][j] * b2[i][j];
			}
		}
		int N = b1.length * (b1.length - 1) / 2;

		double[] br = { Sx1 / N, (Sxx / N - Sx1 / N * Sx2 / N) };
		return br;
	}

	private void HESingle() {
		SimpleRegression sReg = new SimpleRegression();

		int cnt = 0;
		for (int i = 0; i < Y[0].length; i++) {
			for (int j = 0; j < i; j++) {
				double yy = 0;
				if (heBArgs.isSD()) {
					yy = (Y[0][i] - Y[0][j]) * (Y[1][i] - Y[1][j]);
				} else if (heBArgs.isCP()) {
					yy = (Y[0][i] * Y[1][j] + Y[1][i] * Y[0][j]) / 2;
				} else {
					yy = (Y[0][i] + Y[0][j]) * (Y[1][i] + Y[1][j]);
				}
				if (heBArgs.isGRMcut() && A3[0][i][j] > heBArgs.getGRMcutoff())
					continue;
				sReg.addData(A3[0][i][j], yy);
				cnt++;
			}
		}

		Logger.printUserLog("ST " + sReg.getTotalSumSquares() / sReg.getN());
		Logger.printUserLog(cnt + " observations for HE regression.");
		beta = new double[2];
		beta_v = new double[2][2];

		beta[0] = sReg.getIntercept();
		beta[1] = sReg.getSlope();
		beta_v[0][0] = sReg.getInterceptStdErr();
		beta_v[1][1] = sReg.getSlopeStdErr();
	}

	private void HESingleJackknife() {
		double[][] jBeta_ = new double[Y[0].length][2];
		double[] jBeta = new double[2];
		double[] jBetaSe = new double[2];
		for (int jk = 0; jk < Y[0].length; jk++) {
			SimpleRegression sReg = new SimpleRegression();

			for (int i = 0; i < Y[1].length; i++) {
				if (i == jk) {
					continue;
				}
				for (int j = 0; j < i; j++) {
					double yy = 0;
					if (j == jk) {
						continue;
					}
					if (heBArgs.isSD()) {
						yy = (Y[0][i] - Y[0][j]) * (Y[1][i] - Y[1][j]);
					} else if (heBArgs.isCP()) {
						yy = (Y[0][i] * Y[1][j] + Y[1][i] * Y[0][j]) / 2;
					} else {
						yy = (Y[0][i] + Y[0][j]) * (Y[1][i] + Y[1][j]) / 2;
					}
					if (heBArgs.isGRMcut() && A3[0][i][j] > heBArgs.getGRMcutoff())
						continue;
					sReg.addData(A3[0][i][j], yy);
				}
			}
			// Logger.printUserLog(cnt + " observations for HE regression -> Jackknife " +
			// (jk+1) + " " + sReg.getSlope());
			jBeta_[jk][0] = Y[0].length * beta[0] - (Y[0].length - 1) * sReg.getIntercept();
			jBeta_[jk][1] = Y[0].length * beta[1] - (Y[0].length - 1) * sReg.getSlope();
		}

		double s1 = 0, s2 = 0;
		for (int jk = 0; jk < jBeta_.length; jk++) {
			s1 += jBeta_[jk][0];
			s2 += jBeta_[jk][1];
		}
		jBeta[0] = s1 / jBeta_.length;
		jBeta[1] = s2 / jBeta_.length;
		double ss1 = 0, ss2 = 0;
		for (int jk = 0; jk < jBeta_.length; jk++) {
			ss1 += (jBeta_[jk][0] - jBeta[0]) * (jBeta_[jk][0] - jBeta[0]);
			ss2 += (jBeta_[jk][1] - jBeta[1]) * (jBeta_[jk][1] - jBeta[1]);
		}
		jBetaSe[0] = Math.sqrt(ss1 / (jBeta_.length * (jBeta_.length - 1)));
		jBetaSe[1] = Math.sqrt(ss2 / (jBeta_.length * (jBeta_.length - 1)));

		Logger.printUserLog("-------------------------------------------------------------------");
		Logger.printUserLog("Jackknife results: ");
		Logger.printUserLog("JK_Mean\t" + jBeta[0] + "\t" + jBetaSe[0]);
		Logger.printUserLog("JK_Beta1\t" + jBeta[1] + "\t" + jBetaSe[1]);
		Logger.printUserLog("-------------------------------------------------------------------");

	}

	private double[][] readCovar() {
		int[] covIdx = data.getMatchedSubjectIdx(covFileIdx);
		double[][] x = new double[covIdx.length][heBArgs.getCovNumber().length];
		for (int subjectIdx = 0; subjectIdx < covIdx.length; subjectIdx++) {
			x[subjectIdx] = data.getVariable(covFileIdx, covIdx[subjectIdx], heBArgs.getCovNumber());
		}
		return x;
	}

	private double[][] readPhenotypes() {
		int[] pheIdx = data.getMatchedSubjectIdx(pheFileIdx);
		double[][] y = new double[2][pheIdx.length];
		for (int subjectIdx = 0; subjectIdx < y[0].length; subjectIdx++) {
			y[0][subjectIdx] = data.getVariable(pheFileIdx, pheIdx[subjectIdx], heBArgs.getPhenotypeIdx()[0]);
			y[1][subjectIdx] = data.getVariable(pheFileIdx, pheIdx[subjectIdx], heBArgs.getPhenotypeIdx()[1]);
		}
		return y;
	}

	private void adjustY() {
		OLSMultipleLinearRegression mod = new OLSMultipleLinearRegression();
		mod.newSampleData(Y[0], X);
		double[] res = mod.estimateResiduals();
		System.arraycopy(res, 0, Y[0], 0, Y[0].length);

		mod.newSampleData(Y[1], X);
		double[] res1 = mod.estimateResiduals();
		System.arraycopy(res1, 0, Y[1], 0, Y[1].length);
	}

	private double[][][] lineUpGenotype(double[][][] B) {
		int[] subIdx = data.getMatchedSubjectIdx(grmFileIdx);
		double[][][] b = new double[B.length][subIdx.length][subIdx.length];

		for (int i = 0; i < b.length; i++) {
			b[i] = lineUpMatrix(B[i], subIdx);
		}
		return b;
	}

	private double[][] lineUpMatrix(double[][] B, int[] subIdx) {
		double[][] a = new double[subIdx.length][subIdx.length];

		for (int i = 0; i < subIdx.length; i++) {
			for (int j = 0; j < subIdx.length; j++) {
				a[i][j] = a[j][i] = B[subIdx[i]][subIdx[j]];
			}
		}
		return a;
	}

	private double[][][] readGrmList(int numSubjects) {
		double[][][] B = new double[heBArgs.getGrmList().length][numSubjects][numSubjects];
		String[] tokens = null;
		for (int z = 0; z < heBArgs.getGrmList().length; z++) {
			BufferedReader reader = BufferedReader.openGZipFile(heBArgs.getGrmList()[z] + ".grm.gz", "GRM (gzip)");

			for (int i = 0; i < B[z].length; i++) {
				for (int j = 0; j <= i; j++) {
					if ((tokens = reader.readTokens(4)) != null) {
						B[z][i][j] = B[z][j][i] = Double.parseDouble(tokens[3]);
					}
				}
			}
			reader.close();
		}
		return lineUpGenotype(B);
	}

	private void readGrm(int numSubjects) {
		if (heBArgs.getGrmBin() != null) {
			readGrmBin(heBArgs.getGrmBin(), numSubjects);
		} else {
			BufferedReader reader = heBArgs.getGrmText() == null
					? BufferedReader.openGZipFile(heBArgs.getGrmGZ(), "GRM (gzip)")
					: BufferedReader.openTextFile(heBArgs.getGrmText(), "GRM");
			readGrm(reader, numSubjects);
		}
	}

	private void readGrmBin(String fileName, int numSubjects) {
		BinaryInputFile grmBin = new BinaryInputFile(fileName, "GRM (binary)", /* littleEndian */true);
		double[][] B = new double[numSubjects][numSubjects];
		Logger.printUserLog("Constructing A matrix: a " + numSubjects + " X " + numSubjects + " matrix.");
		for (int i = 0; i < B.length; i++) {
			for (int j = 0; j <= i; j++) {
				if (grmBin.available() >= ConstValues.FLOAT_SIZE) {
					B[i][j] = B[j][i] = grmBin.readFloat();
				}
			}
		}
		grmBin.close();
		A3 = new double[1][data.getMatchedSubjectIdx(grmFileIdx).length][data.getMatchedSubjectIdx(grmFileIdx).length];
		A3[0] = lineUpMatrix(B, data.getMatchedSubjectIdx(grmFileIdx));
	}

	private void readGrm(BufferedReader reader, int numSubjects) {
		double[][] B = new double[numSubjects][numSubjects];
		String[] tokens = null;
		for (int i = 0; i < B.length; i++) {
			for (int j = 0; j <= i; j++) {
				if ((tokens = reader.readTokens(4)) != null) {
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

	private double[][][] A3;

	private double[][] Y = null;
	private double[][] X = null;
	private double[] beta = null;
	private double[][] beta_v = null;
	// private double[] jbeta = null;
	private double[] jBeta = null;
	private double[] jBetaSe = null;

	private HERegBCommandArguments heBArgs = null;
	private InputDataSet2 data = null;

	private double[] mA;
	private double[][] vA;

	private DecimalFormat fmt = new DecimalFormat("0.0000");
	private DecimalFormat fmtE = new DecimalFormat("0.00E0");

}
