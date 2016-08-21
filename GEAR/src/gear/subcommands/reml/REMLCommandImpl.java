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
		if (remlArgs.isGRMList())
		{
			InputDataSet data = new InputDataSet();
			data.readSubjectIDFile(remlArgs.getGrmList()[0] + ".grm.id");
			data.readPhenotypeFile(remlArgs.getPhenotypeFile());

			double[] Y = new double[data.getNumberOfSubjects()];
			for(int subjectIdx = 0; subjectIdx < Y.length; subjectIdx++)
			{
				Y[subjectIdx] = data.isPhenotypeMissing(subjectIdx, remlArgs.getPhenotypeIdx()) ? 0 : data.getPhenotype(subjectIdx, remlArgs.getPhenotypeIdx());
			}			

			double[][][] A3 = readGrmList(data.getNumberOfSubjects());
			mlm = new MLM(A3, Y, remlArgs.isMINQUE());
		}
		else
		{
			InputDataSet data = new InputDataSet();
			data.readSubjectIDFile(remlArgs.getGrmID());
			data.readPhenotypeFile(remlArgs.getPhenotypeFile());

			double[] Y = new double[data.getNumberOfSubjects()];
			for(int subjectIdx = 0; subjectIdx < Y.length; subjectIdx++)
			{
				Y[subjectIdx] = data.isPhenotypeMissing(subjectIdx, remlArgs.getPhenotypeIdx()) ? 0 : data.getPhenotype(subjectIdx, remlArgs.getPhenotypeIdx());
			}

			readGrm(data.getNumberOfSubjects());			
			mlm = new MLM(A, Y, remlArgs.isMINQUE());
		}
		mlm.MINQUE();
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
		return B;
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
	REMLCommandArguments remlArgs = null;

}
