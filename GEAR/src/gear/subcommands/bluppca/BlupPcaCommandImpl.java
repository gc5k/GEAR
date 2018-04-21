package gear.subcommands.bluppca;

import gear.ConstValues;
import gear.data.InputDataSet2;
import gear.family.GenoMatrix.GenotypeMatrix;
import gear.family.pedigree.file.SNP;
import gear.family.plink.PLINKParser;
import gear.family.qc.rowqc.SampleFilter;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.BinaryInputFile;
import gear.util.BufferedReader;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.pop.PopStat;

import java.io.PrintStream;
import java.text.DecimalFormat;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;

public class BlupPcaCommandImpl extends CommandImpl {
	private double[][] A;
	private SampleFilter sf;
	private GenotypeMatrix pGM;

	@Override
	public void execute(CommandArguments cmdArgs) {
		BlupPcaCommandArguments blupArgs = (BlupPcaCommandArguments) cmdArgs;

		InputDataSet2 data = new InputDataSet2();
		data.addFile(blupArgs.getGrmID());
		data.addFile(blupArgs.getPhenotypeFile(), blupArgs.getSelectedPhenotype());
		data.LineUpFiles();

		readGrm(blupArgs, data.getNumberOfSubjects());

		PLINKParser pp = PLINKParser.parse(blupArgs);
		sf = new SampleFilter(pp.getPedigreeData());
		pGM = new GenotypeMatrix(sf.getSample(), pp.getMapData());

		double[][] freq = PopStat.calAlleleFrequency(pGM);
		// PopStat.Imputation(gm);

		Logger.printUserLog("Standardizing genotypes...");
		double[][] genoMat = new double[pGM.getNumIndivdial()][pGM.getNumMarker()];
		for (int i = 0; i < genoMat.length; i++) {
			for (int j = 0; j < genoMat[i].length; j++) {
				if (pGM.getAdditiveScore(i, j) == ConstValues.BINARY_MISSING_GENOTYPE) {
					genoMat[i][j] = 0;
				} else {
					genoMat[i][j] = (pGM.getAdditiveScore(i, j) - 2 * freq[j][1])
							/ Math.sqrt(2 * freq[j][1] * (1 - freq[j][1]));
				}
			}
		}

		double[][] blupPC = new double[pGM.getNumMarker()][data.getNumberOfTraits()];

		Logger.printUserLog("Inversing the matrix...");
		RealMatrix grm = new Array2DRowRealMatrix(A);
		RealMatrix grm_Inv = (new LUDecompositionImpl(grm)).getSolver().getInverse();

		// RealMatrix tmp = (new
		// Array2DRowRealMatrix(genoMat)).transpose().multiply(grm_Inv);
		RealMatrix tGenoMat = (new Array2DRowRealMatrix(genoMat)).transpose();

		RealMatrix tmp = new Array2DRowRealMatrix(tGenoMat.getRowDimension(), grm_Inv.getColumnDimension());

		Logger.printUserLog("Revving up the BLUP machine...");

		for (int i = 0; i < tGenoMat.getRowDimension(); i++) {
			for (int j = 0; j < grm_Inv.getColumnDimension(); j++) {
				double f = 0;
				for (int k = 0; k < tGenoMat.getColumnDimension(); k++) {
					if (pGM.getAdditiveScore(k, i) != ConstValues.BINARY_MISSING_GENOTYPE) {
						f += tGenoMat.getEntry(i, k) * grm_Inv.getEntry(k, j);
					}
				}
				tmp.setEntry(i, j, f);
			}
		}

		for (int traitIdx = 0; traitIdx < 1; traitIdx++) {
			Logger.printUserLog("Calculating blup vector[" + (traitIdx + 1) + "].");

			double[] Y = new double[data.getNumberOfSubjects()];
			int[] pIdx = data.getMatchedSubjectIdx(1);
			for (int subjectIdx = 0; subjectIdx < Y.length; subjectIdx++) {
				Y[subjectIdx] = data.getVariable(1, pIdx[subjectIdx], traitIdx);
			}

			RealMatrix B = tmp.multiply(new Array2DRowRealMatrix(Y));

			// Logger.printUserLog("Rescaling the snp effects...");
			for (int j = 0; j < B.getRowDimension(); j++) {
				blupPC[j][traitIdx] = B.getEntry(j, 0) / pGM.getNumMarker();
			}
		}

		PrintStream predictorFile = FileUtil.CreatePrintStream(blupArgs.getOutRoot() + ".blup");

		// Title Line

		predictorFile.print("SNP\tRefAllele");

		for (int i = 0; i < blupArgs.getSelectedPhenotype().length; i++) {
			if (blupArgs.getSelectedPhenotype(i) == 0) {
				predictorFile.print("\tblup_trait" + (i + 1));
			} else {
				predictorFile.print("\tblup_trait" + (blupArgs.getSelectedPhenotype(i) + 1));
			}
		}

		DecimalFormat fmt = new DecimalFormat("#.###E00");

		for (int i = 0; i < pGM.getNumMarker(); i++) {
			SNP snp = pGM.getSNPList().get(i);
			predictorFile.print(snp.getName() + "\t" + snp.getFirstAllele() + "\t");
			for (int j = 0; j < blupPC[i].length; j++) {
				if (j == (blupPC[i].length - 1)) {
					predictorFile.println(fmt.format(blupPC[i][j]));
				} else {
					predictorFile.print(fmt.format(blupPC[i][j]) + "\t");
				}
			}
		}
		predictorFile.close();
	}

	private void readGrm(BlupPcaCommandArguments blupArgs, int numSubjects) {
		if (blupArgs.getGrmBin() != null) {
			readGrmBin(blupArgs.getGrmBin(), numSubjects);
		} else {
			BufferedReader reader = blupArgs.getGrmText() == null
					? BufferedReader.openGZipFile(blupArgs.getGrmGZ(), "GRM (gzip)")
					: BufferedReader.openTextFile(blupArgs.getGrmText(), "GRM");
			readGrm(reader, numSubjects);
		}
	}

	private void readGrmBin(String fileName, int numSubjects) {
		BinaryInputFile grmBin = new BinaryInputFile(fileName, "GRM (binary)", /* littleEndian */true);
		A = new double[numSubjects][numSubjects];
		Logger.printUserLog("Constructing A matrix: a " + numSubjects + " X " + numSubjects + " matrix.");
		for (int i = 0; i < A.length; i++) {
			for (int j = 0; j <= i; j++) {
				if (grmBin.available() >= ConstValues.FLOAT_SIZE) {
					A[i][j] = A[j][i] = grmBin.readFloat();
				}
			}
		}
		grmBin.close();
	}

	private void readGrm(BufferedReader reader, int numSubjects) {
		A = new double[numSubjects][numSubjects];
		String[] tokens = null;
		for (int i = 0; i < A.length; i++) {
			for (int j = 0; j <= i; j++) {
				if ((tokens = reader.readTokens(4)) != null) {
					A[i][j] = A[j][i] = Double.parseDouble(tokens[3]);
				}
			}
		}
		reader.close();
	}

}
