package gear.subcommands.metawatchdog.decrypt;

import java.io.PrintStream;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ChiSquaredDistributionImpl;
import org.apache.commons.math.stat.StatUtils;
import org.apache.commons.math.stat.regression.SimpleRegression;

import gear.ConstValues;
import gear.data.PhenotypeFile;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.BinaryInputFile;
import gear.util.FileUtil;
import gear.util.Logger;

public class MetaWatchdog2CommandImpl extends CommandImpl
{
	@Override
	public void execute(CommandArguments cmdArgs)
	{
		mwArgs = (MetaWatchdog2CommandArguments) cmdArgs;

		initNormalizedScores();

		Logger.printUserLog(normScores[0][0].length + " scores in '" + mwArgs.getDataset1() + "'.");
		Logger.printUserLog(normScores[1][0].length + " scores in '" + mwArgs.getDataset2() + "'.");

		if (normScores[0][0].length != normScores[1][0].length)
		{
			Logger.printUserLog("Warning: the number of scores in two sets are differnt. Only the first " + Math.min(normScores[0][0].length, normScores[1][0].length) + " scores in each file will be used.");
		}


		if (!mwArgs.getRegFlag() && !mwArgs.getChisqFlag())
		{
			if (mwArgs.getEncodeFlag())
			{
				readEncodeFile();
			}			
		}
		else 
		{
			Logger.printUserLog("Encode parameters has been masked manually.");
		}

		if (mwArgs.getChisqFlag())
		{
			Logger.printUserLog("Chi-sq method is used to detect overlapping individuals.");
			Logger.printUserLog("Chi-sq cutoff: " + mwArgs.getChisqCutoff());
			chisqMethod();
		}
		else if (mwArgs.getRegFlag())
		{
			Logger.printUserLog("Regression method is used to detect overlapping individuals.");
			Logger.printUserLog("Regression coefficient cut off: " + mwArgs.getRegB());
			regressionMethod();
		}
	}

	private void chisqMethod()
	{
		PrintStream predictorFile = FileUtil.CreatePrintStream(mwArgs.getOutRoot() + ".mw");

		int numScoreCols = Math.min(normScores[0][0].length, normScores[1][0].length);

		int cntSimilarPairs = 0, cntTotalPairs = 0;
		ChiSquaredDistributionImpl chiDis = new ChiSquaredDistributionImpl(numScoreCols);
		for (int i = 0; i < normScores[0].length; i++)  // # of subjects in data 1
		{
			if (isSubjectIncluded[0][i])
			{
				for (int j = 0; j < normScores[1].length; j++)  // # of subjects in data 2
				{
					if (isSubjectIncluded[1][j])
					{
						double chi = 0;
						for (int k = 0; k < numScoreCols; ++k)
						{
							chi += (normScores[0][i][k] - normScores[1][j][k]) * (normScores[0][i][k] - normScores[1][j][k]) / 2;
						}
						double pchi = 1;
						try
						{
							pchi = chiDis.cumulativeProbability(chi);
						}
						catch (MathException e)
						{
							e.printStackTrace();
						}

						if (chi <= mwArgs.getChisqCutoff())
						{
							String entry = "";
							entry += scores[0].getSubjectID(i).getFamilyID() + "\t" + scores[0].getSubjectID(i).getIndividualID() + "\t";
							entry += scores[1].getSubjectID(j).getFamilyID() + "\t" + scores[1].getSubjectID(j).getIndividualID() + "\t";
							entry += chi + "\t";
							entry += pchi;
							predictorFile.println(entry);
							++cntSimilarPairs;
						}
						++cntTotalPairs;
					}
				}
			}
		}
		predictorFile.close();
		Logger.printUserLog("In total " + cntTotalPairs + " pairs were compared.");
		Logger.printUserLog("In total " + cntSimilarPairs + " similar pairs were detected.");
		if (cntSimilarPairs > 0)
		{
			Logger.printUserLog("Detected pairs have been save in '" + mwArgs.getOutRoot() + ".mw'.");
		}
	}

	private void regressionMethod()
	{
		PrintStream predictorFile = FileUtil.CreatePrintStream(mwArgs.getOutRoot() + ".mw");

		int numScoreCols = Math.min(normScores[0][0].length, normScores[1][0].length); 

		int cntSimilarPairs = 0, cntTotalPairs = 0;
		for (int i = 0; i < normScores[0].length; i++)  // # of subjects in data 1
		{
			if (isSubjectIncluded[0][i])
			{
				for (int j = 0; j < normScores[1].length; j++)  // # of subjects in data 2
				{
					if (isSubjectIncluded[1][j])
					{
						double[][] dat = new double[numScoreCols][2];
						for (int k = 0; k < numScoreCols; ++k)
						{
							dat[k][0] = normScores[0][i][k];
							dat[k][1] = normScores[1][j][k];
						}
						SimpleRegression sr = new SimpleRegression();
						sr.addData(dat);
						double b = sr.getSlope();
		
						if (b > mwArgs.getRegB())
						{
							String entry = "";
							entry += scores[0].getSubjectID(i).getFamilyID() + "\t" + scores[0].getSubjectID(i).getIndividualID() + "\t";
							entry += scores[1].getSubjectID(j).getFamilyID() + "\t" + scores[1].getSubjectID(j).getIndividualID() + "\t";
							entry += b + "\t" + sr.getSlopeStdErr() + "\t" + sr.getN();
							predictorFile.println(entry);
							++cntSimilarPairs;
						}
						++cntTotalPairs;
					}
				}
			}
		}
		predictorFile.close();
		Logger.printUserLog("In total " + cntTotalPairs + " pairs were compared.");
		Logger.printUserLog("In total " + cntSimilarPairs + " similar pairs were detected.");
		if (cntSimilarPairs > 0)
		{
			Logger.printUserLog("Detected pairs have been save in '" + mwArgs.getOutRoot() + ".mw'.");
		}
	}

	private void initNormalizedScores()
	{
		initNormalizedScores(0, mwArgs.getDataset1());
		initNormalizedScores(1, mwArgs.getDataset2());
	}

	private void initNormalizedScores(int dataIdx, String dataFile)
	{
		PhenotypeFile scores = this.scores[dataIdx] = new PhenotypeFile(dataFile, ConstValues.HAS_HEADER);
		double[][] normScores = this.normScores[dataIdx] = new double[scores.getNumberOfSubjects()][scores.getNumberOfTraits()];
		boolean[] isSubjectIncluded = this.isSubjectIncluded[dataIdx] = new boolean[scores.getNumberOfSubjects()];
		PrintStream pstrm = FileUtil.CreatePrintStream(mwArgs.getOutRoot() + ".ignored." + dataFile);
		for (int subjectIdx = 0; subjectIdx < scores.getNumberOfSubjects(); ++subjectIdx)
		{
			isSubjectIncluded[subjectIdx] = true;
			for (int scoreCol = 0; scoreCol < scores.getNumberOfTraits() && isSubjectIncluded[subjectIdx]; ++scoreCol)
			{
				if (scores.isMissing(subjectIdx, scoreCol))
				{
					isSubjectIncluded[subjectIdx] = false;
					pstrm.println(scores.getSubjectID(subjectIdx));
					break;
				}
				normScores[subjectIdx][scoreCol] = scores.getPhenotype(subjectIdx, scoreCol);
			}
			if (isSubjectIncluded[subjectIdx])
			{
				normScores[subjectIdx] = StatUtils.normalize(normScores[subjectIdx]);
			}
		}
		pstrm.close();
	}

	private void readEncodeFile()
	{
		BinaryInputFile file = new BinaryInputFile(mwArgs.getEncodeFile(), "encode");
		double K = file.readDouble();
		long seed = file.readLong();
		double alpha = file.readDouble();
		int tests = file.readInt();

		double beta = file.readDouble();
		double b = file.readDouble();
		double q = file.readDouble();
		int method = file.readInt();
		file.close();

		if (method == 0) // chisq
		{
			Logger.printUserLog("Encode file set the chisq q value to " + q + ".");
			mwArgs.setChisq(q);
		}
		else
		{
			Logger.printUserLog("Encode file set the regression b value to " + b + ".");
			mwArgs.setRegB(b);
		}
	}

	private PhenotypeFile[] scores = new PhenotypeFile[2];
	private double[][][] normScores = new double[2][][];
	private boolean[][] isSubjectIncluded = new boolean[2][];
	private MetaWatchdog2CommandArguments mwArgs;
}
