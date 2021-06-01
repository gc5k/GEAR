package gear.subcommands.grmstat;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;

import gear.ConstValues;
import gear.data.InputDataSet2;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.BinaryInputFile;
import gear.util.BufferedReader;
import gear.util.Logger;
import gear.util.NewIt;

public class GRMStatCommandImpl extends CommandImpl {

	@Override
	public void execute(CommandArguments cmdArgs) {
		gsArgs = (GRMStatCommandArguments) cmdArgs;
		readGrm();
	}

	private void readGrm() {
		if (gsArgs.getGrmBin() != null) {
//			readGrmBin();
			Logger.printUserLog("grm-bin is not supported.");
			System.exit(1);
		} else {
			BufferedReader reader = gsArgs.getGrmText() == null
					? BufferedReader.openGZipFile(gsArgs.getGrmGZ(), "GRM (gzip)")
					: BufferedReader.openTextFile(gsArgs.getGrmText(), "GRM");
			readGrm(reader);
		}
	}

	private void readGrmBin() {
		
		BufferedReader reader = BufferedReader.openTextFile(gsArgs.getGrmID(), "GRM id file");
		String[] tokens = null;
		int sampleCnt= 0;
		while( (tokens = reader.readTokens(2)) != null) {
			sampleCnt++;
		}
		
		reader.close();

		ArrayList<Integer> Diag = NewIt.newArrayList();
		for (int i = 1; i <= sampleCnt; i++) {
			Diag.add(i * (i+1) / 2);
		}

		BinaryInputFile grmBin = new BinaryInputFile(gsArgs.getGrmBin(), "GRM (binary)", /* littleEndian */true);

		double grmMean = 0;
		double grmSq = 0;
		int cnt = 0;

		double grmDMean = 0;
		double grmDSq = 0;
		int cntD = 0;

		int Cnt = 0;
		if (grmBin.available() >= ConstValues.DOUBLE_SIZE) {
			Cnt++;
			double s = grmBin.readDouble();
			if (Collections.binarySearch(Diag, Cnt) >= 0) {
				grmDMean += s;
				grmDSq += s*s;
				cntD++;
			} else {
				grmMean += s;
				grmSq += s*s;
				cnt++;
			}
		}
		grmBin.close();
		
		Logger.printUserLog("Read " + cntD + " diagonal elements.");
		grmDMean /= cntD;
		grmDSq /= cntD;
		double grmDSD = (grmDSq - grmDMean * grmDMean) * cntD / (cntD - 1);

		DecimalFormat df = new DecimalFormat("0.0000");
		DecimalFormat dfE = new DecimalFormat("0.00E0");
		if (Math.abs(grmDMean) > 0.0001) {
			Logger.printUserLog("Mean of the diagonal genetic relatedness is: " + df.format(grmDMean));
		} else {
			Logger.printUserLog("Mean of the diagonal genetic relatedness is: " + dfE.format(grmDMean));
		}

		if (Math.abs(grmDSD) > 0.0001) {
			Logger.printUserLog("Sampling variance of the diagonal genetic relatedness is: " + df.format(grmDSD));
		} else {
			Logger.printUserLog("Sampling variance of the diagonal genetic relatedness is: " + dfE.format(grmDSD));
		}


		Logger.printUserLog("Read " + cnt + " off-diagonal elements.");
		grmMean /= cnt;
		grmSq /= cnt;
		double Effective_sample = -1 / grmMean + 1;
		double grmSD = (grmSq - grmMean * grmMean) * cnt / (cnt - 1);
		double Effective_marker = 1 / grmSD;

		if (Math.abs(grmMean) > 0.0001) {
			Logger.printUserLog("Mean of genetic relatedness is: " + df.format(grmMean));
		} else {
			Logger.printUserLog("Mean of genetic relatedness is: " + dfE.format(grmMean));
		}

		if (Math.abs(grmSD) > 0.0001) {
			Logger.printUserLog("Sampling variance of genetic relatedness is: " + df.format(grmSD));
		} else {
			Logger.printUserLog("Sampling variance of genetic relatedness is: " + dfE.format(grmSD));
		}

		if (Math.abs(Effective_marker) > 0.0001) {
			Logger.printUserLog("Effective sample size is: " + df.format(Effective_sample));
		} else {
			Logger.printUserLog("Effective sample size is: " + dfE.format(Effective_sample));
		}
		
		if (Math.abs(Effective_marker) > 0.0001) {
			Logger.printUserLog("Effective number of genome segments is: " + df.format(Effective_marker));
		} else {
			Logger.printUserLog("Effective number of genome segments is: " + dfE.format(Effective_marker));
		}

	}

	private void readGrm(BufferedReader reader) {
		String[] tokens = null;
		double grmMean = 0;
		double grmSq = 0;
		int cnt = 0;

		double grmDMean = 0;
		double grmDSq = 0;
		int cntD = 0;

		while((tokens = reader.readTokens(4)) != null) {
			int id1 = Integer.parseInt(tokens[0]);
			int id2 = Integer.parseInt(tokens[1]);
			double s = Double.parseDouble(tokens[3]);

			if (id1 == id2) {
				grmDMean += s;
				grmDSq += s*s;
				cntD++;
			} else {
				grmMean += s;
				grmSq += s*s;
				cnt++;
			}
		}
		reader.close();

		Logger.printUserLog("Read " + cntD + " diagonal elements.");
		grmDMean /= cntD;
		grmDSq /= cntD;
		double grmDSD = (grmDSq - grmDMean * grmDMean) * cntD / (cntD - 1);

		DecimalFormat df = new DecimalFormat("0.0000");
		DecimalFormat dfE = new DecimalFormat("0.00E0");
		if (Math.abs(grmDMean) > 0.0001) {
			Logger.printUserLog("Mean of the diagonal genetic relatedness is: " + df.format(grmDMean));
		} else {
			Logger.printUserLog("Mean of the diagonal genetic relatedness is: " + dfE.format(grmDMean));
		}

		if (Math.abs(grmDSD) > 0.0001) {
			Logger.printUserLog("Sampling variance of the diagonal genetic relatedness is: " + df.format(grmDSD));
		} else {
			Logger.printUserLog("Sampling variance of the diagonal genetic relatedness is: " + dfE.format(grmDSD));
		}


		Logger.printUserLog("Read " + cnt + " off-diagonal elements.");
		grmMean /= cnt;
		grmSq /= cnt;
		double Effective_sample = -1 / grmMean + 1;
		double grmSD = (grmSq - grmMean * grmMean) * cnt / (cnt - 1);
		double Effective_marker = 1 / grmSD;

		if (Math.abs(grmMean) > 0.0001) {
			Logger.printUserLog("Mean of genetic relatedness is: " + df.format(grmMean));
		} else {
			Logger.printUserLog("Mean of genetic relatedness is: " + dfE.format(grmMean));
		}

		if (Math.abs(grmSD) > 0.0001) {
			Logger.printUserLog("Sampling variance of genetic relatedness is: " + df.format(grmSD));
		} else {
			Logger.printUserLog("Sampling variance of genetic relatedness is: " + dfE.format(grmSD));
		}

		if (Math.abs(Effective_marker) > 0.0001) {
			Logger.printUserLog("Effective sample size is: " + df.format(Effective_sample));
		} else {
			Logger.printUserLog("Effective sample size is: " + dfE.format(Effective_sample));
		}
		
		if (Math.abs(Effective_marker) > 0.0001) {
			Logger.printUserLog("Effective number of genome segments is: " + df.format(Effective_marker));
		} else {
			Logger.printUserLog("Effective number of genome segments is: " + dfE.format(Effective_marker));
		}
	}

	private double[][] readGrmGZMatrix() {
		InputDataSet2 data = new InputDataSet2();
		data.addFile(gsArgs.getGrmID());
		int numSubjects = data.getFileSampleSize(grmFileIdx);
		double[][] B = new double[numSubjects][numSubjects];
		String[] tokens = null;
		BufferedReader reader = BufferedReader.openGZipFile(gsArgs.getGrmGZ(), "GRM (gzip)");

		for (int i = 0; i < B.length; i++) {
			for (int j = 0; j <= i; j++) {
				if ((tokens = reader.readTokens(4)) != null) {
					B[i][j] = B[j][i] = Double.parseDouble(tokens[3]);
				}
			}
		}
		reader.close();
		return B;
	}

	private GRMStatCommandArguments gsArgs = null;
	private int grmFileIdx = 0;

}
