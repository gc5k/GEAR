package gear.subcommands.metawatchdog.decrypt;

import java.io.PrintStream;

import org.apache.commons.math.stat.StatUtils;
import org.apache.commons.math.stat.regression.SimpleRegression;

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
		MetaWatchdog2CommandArguments mwArgs = (MetaWatchdog2CommandArguments)cmdArgs;
		
		Logger.printUserLog("Cutoff: " + mwArgs.getCutoff());
		
		PhenotypeFile phe1 = new PhenotypeFile(mwArgs.getDataset1(), /*hasHeaders*/true);
		PhenotypeFile phe2 = new PhenotypeFile(mwArgs.getDataset2(), /*hasHeaders*/true);
		
		PrintStream predictorFile = FileUtil.CreatePrintStream(mwArgs.getOutRoot() + ".watchdog");

		int test = Math.min(phe1.getNumberOfTraits(), phe2.getNumberOfTraits()); 

		int cnt = 0;
		for (int i = 0; i < phe1.getNumberOfSubjects(); i++)
		{
			for (int j = 0; j < phe2.getNumberOfSubjects(); j++)
			{
				// TODO: Move the normalization out of the for-loops to speed-up.
				double[][] normPhes = new double[2][test];
				for (int k = 0; k < test; ++k)
				{
					normPhes[0][k] = phe1.getPhenotype(i, k);
					normPhes[1][k] = phe2.getPhenotype(j, k);
				}
				normPhes[0] = StatUtils.normalize(normPhes[0]);
				normPhes[1] = StatUtils.normalize(normPhes[1]);
				
				double[][] dat = new double[test][2];
				for (int k = 0; k < test; ++k)
				{
					dat[k][0] = normPhes[0][k];
					dat[k][1] = normPhes[1][k];
				}
				SimpleRegression sr = new SimpleRegression();
				sr.addData(dat);
				double b = sr.getSlope();

				if (b > mwArgs.getCutoff())
				{
					String entry = "";
					entry += phe1.getSubjectID(i).getFamilyID() + "\t" + phe1.getSubjectID(i).getIndividualID() + "\t";
					entry += phe2.getSubjectID(j).getFamilyID() + "\t" + phe2.getSubjectID(j).getIndividualID() + "\t";
					entry += b + "\t" + sr.getSlopeStdErr() + "\t" + sr.getN();
					predictorFile.println(entry);
					cnt++;
				}
			}
		}
		predictorFile.close();
		Logger.printUserLog("In total " + phe1.getNumberOfSubjects() * phe2.getNumberOfSubjects() + " pairs were compared.");
		Logger.printUserLog("In total " + cnt + " similar pairs were detected.");
	}
}
