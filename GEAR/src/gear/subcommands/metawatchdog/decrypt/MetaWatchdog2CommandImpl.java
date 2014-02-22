package gear.subcommands.metawatchdog.decrypt;

import java.io.PrintStream;

import org.apache.commons.math.stat.StatUtils;
import org.apache.commons.math.stat.regression.SimpleRegression;

import gear.ConstValues;
import gear.data.PhenotypeFile;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.FileUtil;
import gear.util.Logger;

public class MetaWatchdog2CommandImpl extends CommandImpl
{
	@Override
	public void execute(CommandArguments cmdArgs)
	{
		mwArgs = (MetaWatchdog2CommandArguments)cmdArgs;
		
		Logger.printUserLog("Cutoff: " + mwArgs.getCutoff());
		
		initNormalizedScores();
		
		PrintStream predictorFile = FileUtil.CreatePrintStream(mwArgs.getOutRoot() + ".watchdog");

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
		
						if (b > mwArgs.getCutoff())
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
	
	private PhenotypeFile[] scores = new PhenotypeFile[2];
	private double[][][] normScores = new double[2][][];
	private boolean[][] isSubjectIncluded = new boolean[2][];
	MetaWatchdog2CommandArguments mwArgs;
}
