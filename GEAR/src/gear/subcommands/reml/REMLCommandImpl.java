package gear.subcommands.reml;

import gear.ConstValues;
import gear.data.InputDataSet;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.BinaryInputFile;
import gear.util.BufferedReader;
import gear.util.Logger;
import gear.util.stat.MLM;

public class REMLCommandImpl extends CommandImpl {

	@Override
	public void execute(CommandArguments cmdArgs) 
	{
		remlArgs = (REMLCommandArguments)cmdArgs;

		MLM mlm = null;
		double[] Y = null;
		double[][] X = null;
		if (remlArgs.isGRMList())
		{
			if (remlArgs.getCovFile() == null)
			{
				data = new InputDataSet(remlArgs.getGrmList()[0] + ".grm.id", remlArgs.getPhenotypeFile(), remlArgs.getPhenotypeIdx());
			}
			else
			{
				data = new InputDataSet(remlArgs.getGrmList()[0] + ".grm.id", remlArgs.getPhenotypeFile(), remlArgs.getCovFile(), remlArgs.getPhenotypeIdx(), remlArgs.getCovNumber());
				X = readCovar();
			}

			Y = readPhenotype();
			double[][][] A3 = readGrmList(data.getSubjectFileSampleSize());
			Logger.printUserLog(A3.length + " variance components included.");
			Logger.printUserLog("");

			if (remlArgs.getCovFile() == null)
			{
				mlm = new MLM(A3, Y, remlArgs.isMINQUE());				
			}
			else
			{
				mlm = new MLM(A3, X, Y, remlArgs.isMINQUE(), remlArgs.getCovNumber());				
			}
		}
		else
		{
			if (remlArgs.getCovFile() == null)
			{
				data = new InputDataSet(remlArgs.getGrmID(), remlArgs.getPhenotypeFile(), remlArgs.getPhenotypeIdx());				
			}
			else
			{
				data = new InputDataSet(remlArgs.getGrmID(), remlArgs.getPhenotypeFile(), remlArgs.getCovFile(), remlArgs.getPhenotypeIdx(), remlArgs.getCovNumber());				
				X = readCovar();
			}

			Y = readPhenotype();

			readGrm(data.getSubjectFileSampleSize());
			Logger.printUserLog("1 variance component included.");
			Logger.printUserLog("");

			if (remlArgs.getCovFile() == null)
			{
				mlm = new MLM(A, Y, remlArgs.isMINQUE());				
			}
			else
			{
				mlm = new MLM(A, X, Y, remlArgs.isMINQUE(), remlArgs.getCovNumber());
			}
		}

		mlm.MINQUE();
		mlm.printVC();
	}

	private double[][] readCovar()
	{
		int[] covIdx = data.getMatchedCovSubIdx();
		double[][] x = new double[covIdx.length][remlArgs.getCovNumber().length];
		for (int subjectIdx = 0; subjectIdx < covIdx.length; subjectIdx++)
		{
			for (int j = 0; j < remlArgs.getCovNumber().length; j++)
			{
				x[subjectIdx][j] = data.getCovariate(subjectIdx, remlArgs.getCovNumber()[j]);
			}
		}
		return x;
	}

	private double[] readPhenotype()
	{
		int[] pheIdx = data.getMatchedPheSubIdx();
		double[] y = new double[pheIdx.length];
		for(int subjectIdx = 0; subjectIdx < y.length; subjectIdx++)
		{
			y[subjectIdx] = data.getPhenotype(pheIdx[subjectIdx], remlArgs.getPhenotypeIdx());
		}
		return y;
	}
	
	private double[][][] lineUpGenotype(double[][][] B)
	{
		int[] subIdx = data.getMatchedSubIdx();
		double[][][] b = new double[B.length][subIdx.length][subIdx.length];
		
		for (int i = 0; i < b.length; i++)
		{
			for (int j = 0; j < subIdx.length; j++)
			{
				for (int k = 0; k <= j; k++)
				{
					b[i][j][k] = b[i][k][j] = B[i][subIdx[j]][subIdx[k]];
				}
			}
		}
		return b;
	}

	private void lineUpGenotype(double[][] B)
	{
		int[] subIdx = data.getMatchedSubIdx();
		A = new double[subIdx.length][subIdx.length];
		
		for (int i = 0; i < subIdx.length; i++)
		{
			for (int j = 0; j < subIdx.length; j++)
			{
				A[i][j] = A[j][i] = B[subIdx[i]][subIdx[j]];
			}
		}
	}

	private double[][][] readGrmList(int numSubjects)
	{
		double[][][] B = new double[remlArgs.getGrmList().length][numSubjects][numSubjects];
		String[] tokens = null;
		for (int z = 0; z < remlArgs.getGrmList().length; z++)
		{
			BufferedReader reader = BufferedReader.openGZipFile(remlArgs.getGrmList()[z] + ".grm.gz", "GRM (gzip)");

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
		if (remlArgs.getGrmBin() != null)
		{
			readGrmBin(remlArgs.getGrmBin(), numSubjects);
		}
		else
		{
			BufferedReader reader = remlArgs.getGrmText() == null ?
					BufferedReader.openGZipFile(remlArgs.getGrmGZ(), "GRM (gzip)") :
					BufferedReader.openTextFile(remlArgs.getGrmText(), "GRM");
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
		lineUpGenotype(B);
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
		lineUpGenotype(B);
	}

	private double[][] A;
	private REMLCommandArguments remlArgs = null;
	private InputDataSet data = null;
}
