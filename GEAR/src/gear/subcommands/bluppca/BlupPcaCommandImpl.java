package gear.subcommands.bluppca;

import gear.ConstValues;
import gear.data.InputDataSet;
import gear.family.pedigree.file.MapFile;
import gear.family.pedigree.file.SNP;
import gear.family.plink.PLINKParser;
import gear.family.popstat.GenotypeMatrix;
import gear.family.qc.rowqc.SampleFilter;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.sumstat.qc.rowqc.SumStatQC;
import gear.util.BinaryInputFile;
import gear.util.BufferedReader;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.pop.PopStat;

import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;

public class BlupPcaCommandImpl extends CommandImpl
{
	@Override
	public void execute(CommandArguments cmdArgs)
	{
		BlupPcaCommandArguments blupArgs = (BlupPcaCommandArguments)cmdArgs;

		InputDataSet data = new InputDataSet();
		data.readSubjectIDFile(blupArgs.getGrmID());
		data.readPhenotypeFile(blupArgs.getPhenotypeFile());
		
		readGrm(blupArgs, data.getNumberOfSubjects());

		PLINKParser pp = PLINKParser.parse(blupArgs);
		sf = new SampleFilter(pp.getPedigreeData(), pp.getMapData());
		ssQC = new SumStatQC(pp.getPedigreeData(), pp.getMapData(), sf);
		mapFile = ssQC.getMapFile();
		gm = new GenotypeMatrix(ssQC.getSample());

		double[][] freq=PopStat.calAlleleFrequency(gm, gm.getNumMarker());
//		PopStat.Imputation(gm);

		Logger.printUserLog("Standardizing genotypes...");
		double[][] genoMat = new double[gm.getNumIndivdial()][gm.getNumMarker()];
		for(int i = 0; i < genoMat.length; i++)
		{
			for(int j = 0; j < genoMat[i].length; j++)
			{
				if (gm.getAdditiveScore(i, j) == ConstValues.BINARY_MISSING_GENOTYPE)
				{
					genoMat[i][j] = 0;
				}
				else
				{
					if (freq[j][1] < 0.001)
					{
						genoMat[i][j] = 0;	
					}
					else
					{
						genoMat[i][j] = (gm.getAdditiveScore(i, j) - 2 * freq[j][1])/Math.sqrt(2*freq[j][1] * (1-freq[j][1]));						
					}
				}
			}
		}

		double[][] blupPC = new double[gm.getNumMarker()][data.getNumberOfTraits()];

		Logger.printUserLog("Inversing the matrix...");
		RealMatrix grm = new Array2DRowRealMatrix(A);
		RealMatrix grm_Inv = (new LUDecompositionImpl(grm)).getSolver().getInverse();

//		RealMatrix tmp = (new Array2DRowRealMatrix(genoMat)).transpose().multiply(grm_Inv);
		RealMatrix tGenoMat = (new Array2DRowRealMatrix(genoMat)).transpose();

		RealMatrix tmp = new Array2DRowRealMatrix(tGenoMat.getRowDimension(), grm_Inv.getColumnDimension());

		Logger.printUserLog("Revving up the BLUP machine...");

		for(int i = 0; i < tGenoMat.getRowDimension(); i++)
		{
			for(int j = 0; j < grm_Inv.getColumnDimension(); j++)
			{
				double f = 0;
				for(int k = 0; k < tGenoMat.getColumnDimension(); k++)
				{
					if (gm.getAdditiveScore(k, i) != ConstValues.BINARY_MISSING_GENOTYPE)
					{
						f += tGenoMat.getEntry(i, k) * grm_Inv.getEntry(k, j);
					}
				}
				tmp.setEntry(i, j, f);
			}
		}

		for(int traitIdx = 0; traitIdx < data.getNumberOfTraits(); traitIdx++)
		{
			Logger.printUserLog("Calculating blup vector[" + (traitIdx+1) + "].");

			double[] Y = new double[data.getNumberOfSubjects()];
			for(int subjectIdx = 0; subjectIdx < Y.length; subjectIdx++)
			{
				Y[subjectIdx] = data.isPhenotypeMissing(subjectIdx, traitIdx) ? 0 : data.getPhenotype(subjectIdx, traitIdx);
			}
			RealMatrix B = tmp.multiply(new Array2DRowRealMatrix(Y));

//			Logger.printUserLog("Rescaling the snp effects...");
			for(int j = 0; j < B.getRowDimension(); j++)
			{
				blupPC[j][traitIdx] = B.getEntry(j, 0) / gm.getNumMarker();
			}
		}

		PrintStream predictorFile = FileUtil.CreatePrintStream(blupArgs.getOutRoot() + ".blup");

		// Title Line
		ArrayList<SNP> snpList = mapFile.getMarkerList();

		predictorFile.print("SNP\tRefAllele");
		for(int traitIdx = 0; traitIdx < data.getNumberOfTraits(); traitIdx++)
		{
			if (traitIdx == (data.getNumberOfTraits() - 1)) 
			{
				predictorFile.println("\tBLUP" + (traitIdx+1));				
			}
			else
			{
				predictorFile.print("\tBLUP" + (traitIdx+1));
			}
		}

		DecimalFormat fmt = new DecimalFormat("#.###E00");

		for(int i = 0; i < gm.getNumMarker(); i++)
		{
			SNP snp = snpList.get(i);
			predictorFile.print(snp.getName() + "\t" + snp.getFirstAllele() + "\t");
			for(int j = 0; j < blupPC[i].length; j++)
			{
				if (j == (blupPC[i].length - 1))
				{
					predictorFile.println(fmt.format(blupPC[i][j]));
				}
				else
				{
					predictorFile.print(fmt.format(blupPC[i][j])+"\t");
				}
			}
		}
		predictorFile.close();
	}

	private void readGrm(BlupPcaCommandArguments blupArgs, int numSubjects)
	{
		if (blupArgs.getGrmBin() != null)
		{
			readGrmBin(blupArgs.getGrmBin(), numSubjects);
		}
		else
		{
			BufferedReader reader = blupArgs.getGrmText() == null ?
					BufferedReader.openGZipFile(blupArgs.getGrmGZ(), "GRM (gzip)") :
					BufferedReader.openTextFile(blupArgs.getGrmText(), "GRM");
			readGrm(reader, numSubjects);
		}
	}

	private void readGrmBin(String fileName, int numSubjects)
	{
		BinaryInputFile grmBin = new BinaryInputFile(fileName, "GRM (binary)", /*littleEndian*/true);
		A = new double[numSubjects][numSubjects];
		Logger.printUserLog("Constructing A matrix: a " + numSubjects + " X " + numSubjects + " matrix.");
		for (int i = 0; i < A.length; i++) 
		{
			for (int j = 0; j <= i; j++)
			{
				if (grmBin.available() >= ConstValues.FLOAT_SIZE)
				{
					A[i][j] = A[j][i] = grmBin.readFloat();
				}
			}
		}
		grmBin.close();
	}

	private void readGrm(BufferedReader reader, int numSubjects)
	{
		A = new double[numSubjects][numSubjects];
		String[] tokens = null;
		for (int i = 0; i < A.length; i++)
		{
			for (int j = 0; j <= i; j++)
			{
				if ((tokens = reader.readTokens(4)) != null) 
				{
					A[i][j] = A[j][i] = Double.parseDouble(tokens[3]);
				}
			}
		}
		reader.close();
	}

	private double[][] A;

	private MapFile mapFile;
	private SampleFilter sf;
	private SumStatQC ssQC;
	private GenotypeMatrix gm;
}
