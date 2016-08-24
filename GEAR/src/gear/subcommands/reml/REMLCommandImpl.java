package gear.subcommands.reml;

import gear.ConstValues;
import gear.data.InputDataSet2;
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
		data = new InputDataSet2();
		int fidx = 1;

		if (remlArgs.isGRMList())
		{
			data.addFile(remlArgs.getGrmList()[0]+".grm.id");
			data.addFile(remlArgs.getPhenotypeFile(), remlArgs.getPhenotypeIdx());

			if (remlArgs.getCovFile() != null)
			{
				data.addFile(remlArgs.getCovFile(), remlArgs.getCovNumber());
				fidx++;
				covFileIdx = fidx;
			}
			if (remlArgs.getKeepFile() != null)
			{
				data.addFile(remlArgs.getKeepFile());
				fidx++;
				keepFileIdx = fidx;
			}
			data.LineUpFiles();

			Y = readPhenotype();
			double[][][] A3 = readGrmList(data.getFileSampleSize(grmFileIdx));
			Logger.printUserLog(A3.length + " variance components included.");
			Logger.printUserLog("");

			if (remlArgs.getCovFile() == null)
			{
				mlm = new MLM(A3, Y, remlArgs.isMINQUE());
			}
			else
			{
				X = readCovar();
				mlm = new MLM(A3, X, Y, remlArgs.isMINQUE(), remlArgs.getCovNumber());				
			}
		}
		else
		{
			data.addFile(remlArgs.getGrmID());
			data.addFile(remlArgs.getPhenotypeFile(), remlArgs.getPhenotypeIdx());

			if (remlArgs.getCovFile() != null)
			{
				data.addFile(remlArgs.getCovFile(), remlArgs.getCovNumber());
				fidx++;
				covFileIdx = fidx;
			}
			if (remlArgs.getKeepFile() != null)
			{
				data.addFile(remlArgs.getKeepFile());
				fidx++;
				keepFileIdx = fidx;
			}
			data.LineUpFiles();

			Y = readPhenotype();

			readGrm(data.getFileSampleSize(grmFileIdx));
			Logger.printUserLog("1 variance component included.");
			Logger.printUserLog("");

			if (remlArgs.getCovFile() == null)
			{
				mlm = new MLM(A, Y, remlArgs.isMINQUE());				
			}
			else
			{
				X = readCovar();
				mlm = new MLM(A, X, Y, remlArgs.isMINQUE(), remlArgs.getCovNumber());
			}
		}

		mlm.MINQUE();
		mlm.printVC();
	}

	private double[][] readCovar()
	{
		int[] covIdx = data.getMatchedSubjectIdx(covFileIdx);
		double[][] x = new double[covIdx.length][remlArgs.getCovNumber().length];
		for (int subjectIdx = 0; subjectIdx < covIdx.length; subjectIdx++)
		{
			x[subjectIdx] = data.getVariable(covFileIdx, covIdx[subjectIdx], remlArgs.getCovNumber());
		}
		return x;
	}

	private double[] readPhenotype()
	{
		int[] pheIdx = data.getMatchedSubjectIdx(pheFileIdx);
		double[] y = new double[pheIdx.length];
		for(int subjectIdx = 0; subjectIdx < y.length; subjectIdx++)
		{
			y[subjectIdx] = data.getVariable(pheFileIdx, pheIdx[subjectIdx], remlArgs.getPhenotypeIdx()[0]);
		}
		return y;
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
		A = lineUpMatrix(B, data.getMatchedSubjectIdx(grmFileIdx));
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
		A = lineUpMatrix(B, data.getMatchedSubjectIdx(grmFileIdx));
	}

	private int grmFileIdx = 0;
	private int pheFileIdx = 1;
	private int covFileIdx = 2;
	private int keepFileIdx = 3;

	private double[][] A;
	private REMLCommandArguments remlArgs = null;
	private InputDataSet2 data = null;
}
