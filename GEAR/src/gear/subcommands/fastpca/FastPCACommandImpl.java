package gear.subcommands.fastpca;

import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.EigenDecomposition;
import org.apache.commons.math.linear.EigenDecompositionImpl;
import org.apache.commons.math.linear.QRDecomposition;
import org.apache.commons.math.linear.QRDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;

import gear.data.InputDataSet2;
import gear.data.SubjectID;
import gear.family.GenoMatrix.GenotypeMatrix;
import gear.family.plink.PLINKParser;
import gear.family.qc.rowqc.SampleFilter;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.subcommands.grmA.GRMACommandArguments;
import gear.subcommands.grmA.GRMACommandImpl;
import gear.util.BufferedReader;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.Sample;
import gear.util.stat.PCA.NCGStatUtils;

public class FastPCACommandImpl extends CommandImpl {
	private double[][] A;
	private int SampleSize;
	private SampleFilter sf;
	private GenotypeMatrix pGM;

	// ArrayList<ArrayList<String>> famID = NewIt.newArrayList();
	private FastPCACommandArguments fpcaArgs;
	private InputDataSet2 data;

	double[] EigenValue;
	double[][] EigenVector;

	@Override
	public void execute(CommandArguments cmdArgs) {
		fpcaArgs = (FastPCACommandArguments) cmdArgs;

		GRMACommandArguments grmArgs = new GRMACommandArguments();
		grmArgs.setBFile(fpcaArgs.getBFile());
		grmArgs.setOutRoot(fpcaArgs.getOutRoot());
		grmArgs.setGZ();
		if (fpcaArgs.isAdjVar())
			grmArgs.setAdjVar();

		GRMACommandImpl grmImpl = new GRMACommandImpl();
		grmImpl.execute(grmArgs);

		data = new InputDataSet2();
		data.addFile(fpcaArgs.getOutRoot() + ".grm.id");
		if (fpcaArgs.isKeepFile())
			data.addFile(fpcaArgs.getKeepFile());
		data.LineUpFiles();

		readGrm();

		PLINKParser pp = PLINKParser.parse(fpcaArgs);
		sf = new SampleFilter(pp.getPedigreeData(), cmdArgs);
		pGM = new GenotypeMatrix(sf.getSample(), pp.getMapData(), cmdArgs);

		SampleSize = (int) Math.ceil(A.length * fpcaArgs.getProp());
		Logger.printUserLog("Sampled " + SampleSize + " individuals.");

		// EigenAnalysis();
		FastPCA();
		WriteEigen();
	}

	private void FastPCA() {
		double[][] genoMat = new double[pGM.getNumIndivdial()][pGM.getNumMarker()];
		for (int i = 0; i < genoMat.length; i++) {
			for (int j = 0; j < genoMat[i].length; j++) {
				genoMat[i][j] = pGM.getAdditiveScore(i, j);
			}
		}

		genoMat = NCGStatUtils.standardize(genoMat, true);

		double[][] genoR = new double[SampleSize][pGM.getNumMarker()];
		Sample.setSeed(fpcaArgs.getSeed());
		int[] idx = Sample.SampleIndex(0, pGM.getNumIndivdial() - 1, SampleSize);
		Arrays.sort(idx);

		for (int i = 0; i < idx.length; i++)
			System.arraycopy(genoMat[idx[i]], 0, genoR[i], 0, genoR[i].length);

		genoR = NCGStatUtils.transpose(genoR);

		Array2DRowRealMatrix gmat = new Array2DRowRealMatrix(genoMat);
		Array2DRowRealMatrix gr = new Array2DRowRealMatrix(genoR);
		double[][] genoY1 = gmat.multiply(gr).getData();
		genoY1 = NCGStatUtils.standardize(genoY1, true);

		double[][] genoY = ((new Array2DRowRealMatrix(A)).multiply(new Array2DRowRealMatrix(genoY1))).getData();
		genoY = NCGStatUtils.standardize(genoY, true);

		QRDecomposition qr = new QRDecompositionImpl(new Array2DRowRealMatrix(genoY));
		RealMatrix Q = qr.getQ().getSubMatrix(0, genoY.length - 1, 0, genoY[0].length - 1);

		RealMatrix B = Q.transpose().multiply(new Array2DRowRealMatrix(genoMat));
		RealMatrix S = B.multiply(B.transpose());

		EigenDecomposition ed = new EigenDecompositionImpl(S, 1e-6);
		EigenValue = ed.getRealEigenvalues();
		for (int i = 0; i < EigenValue.length; i++) {
			EigenValue[i] = Math.sqrt(EigenValue[i] / (pGM.getNumIndivdial() - 1));
		}

		EigenVector = Q.multiply(ed.getV()).getData();
	}

	private void readGrm() {
		BufferedReader reader = BufferedReader.openGZipFile(fpcaArgs.getOutRoot() + ".grm.gz", "GRM (gzip)");
		readGrm(reader);
	}

	private void readGrm(BufferedReader reader) {
		double[][] grmMat = new double[data.getFileSampleSize(0)][data.getFileSampleSize(0)];
		String[] tokens = null;
		for (int i = 0; i < grmMat.length; i++) {
			for (int j = 0; j <= i; j++) {
				if ((tokens = reader.readTokens(4)) != null) {
					grmMat[i][j] = grmMat[j][i] = Double.parseDouble(tokens[3]);
				}
			}
		}
		reader.close();

		A = lineUpMatrix(grmMat, data.getMatchedSubjectIdx(0));
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

	private void WriteEigen() {
		DecimalFormat fmt = new DecimalFormat("0.0000");
		DecimalFormat fmtp = new DecimalFormat("0.000E000");

		PrintStream evaWriter = FileUtil.CreatePrintStream(new String(fpcaArgs.getOutRoot() + ".eigenval"));

		for (int i = 0; i < EigenValue.length; i++) {
			if (Math.abs(EigenValue[i]) >= 0.0001) {
				evaWriter.println(fmt.format(EigenValue[i]));
			} else {
				evaWriter.println(fmtp.format(EigenValue[i]));
			}

			if ((i + 1) <= fpcaArgs.getEV()) {
				Logger.printUserLog("The " + (i + 1) + "th eigenvalue is " + EigenValue[i] + ".");
			}
		}
		evaWriter.close();

		PrintStream eveWriter = FileUtil.CreatePrintStream(new String(fpcaArgs.getOutRoot() + ".eigenvec"));

		ArrayList<SubjectID> SID = data.getMatchedSubjectID(0);

		for (int i = 0; i < EigenVector.length; i++) {
			eveWriter.print(SID.get(i).toString() + "\t");
			for (int j = 0; j < fpcaArgs.getEV(); j++) {
				if (Math.abs(EigenVector[i][j]) >= 0.0001) {
					eveWriter.print(fmt.format(EigenVector[i][j]) + "\t");
				} else {
					eveWriter.print(fmtp.format(EigenVector[i][j]) + "\t");
				}
			}
			eveWriter.println();
		}
		eveWriter.close();
	}
}
