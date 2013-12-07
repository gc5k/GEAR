package gear.metawatchdog.decrypt;

import gear.CmdArgs;
import gear.data.InputDataSet;
import gear.data.SubjectID;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;
import gear.util.stat.PCA.NCGStatUtils;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;
import org.apache.commons.math.distribution.ChiSquaredDistribution;
import org.apache.commons.math.distribution.ChiSquaredDistributionImpl;

public class MetaWatchdog
{
	private final String delim = "\\s+";
	private double[][] phe1;
	private ArrayList<String> ID1 = NewIt.newArrayList();
	private String f1 = null;
	private double[][] phe2;
	private ArrayList<String> ID2 = NewIt.newArrayList();
	private String f2 = null;
	private double alpha = 0.05;
	private NormalDistribution nd = new NormalDistributionImpl();
	private NormalDistribution ndp = new NormalDistributionImpl(0, Math.sqrt(2));
	private ChiSquaredDistribution chip;

	private double T = 0;

	private double[][] pheCom;

	public MetaWatchdog()
	{
		f1 = CmdArgs.INSTANCE.set1_file;
		f2 = CmdArgs.INSTANCE.set2_file;
		alpha = CmdArgs.INSTANCE.dog_alpha;
		try
		{
			T = nd.inverseCumulativeProbability(alpha / 2 + 0.5) * Math.sqrt(2);
		}
		catch (MathException e)
		{
			e.printStackTrace();
		}
		Logger.printUserLog("meta-watchdog sequential test.");
		Logger.printUserLog("set1: " + f1);
		Logger.printUserLog("Set2: " + f2);
		Logger.printUserLog("alpha: " + alpha);
		Logger.printUserLog("cutoff: " + T);
	}

	public void Bark()
	{
		readPhenotypes1(f1);
		readPhenotypes2(f2);

		pheCom = new double[phe1.length + phe2.length][phe1[0].length];
		for (int i = 0; i < phe1.length; i++)
		{
			System.arraycopy(phe1[i], 0, pheCom[i], 0, phe1[i].length);
		}

		for (int i = 0; i < phe2.length; i++)
		{
			System.arraycopy(phe2[i], 0, pheCom[phe1.length+i], 0, phe2[i].length);
		}
		pheCom = NCGStatUtils.standardize(pheCom, false);

		for (int i = 0; i < phe1.length; i++)
		{
			System.arraycopy(pheCom[i], 0, phe1[i], 0, phe1[i].length);
		}

		for (int i = 0; i < phe2.length; i++)
		{
			System.arraycopy(pheCom[phe1.length+i], 0, phe2[i], 0, phe2[i].length);
		}

		PrintStream predictorFile = FileUtil
				.CreatePrintStream(CmdArgs.INSTANCE.out + ".watchdog.seq");
		int test = phe1[0].length <= phe2[0].length ? phe1[0].length
				: phe2[0].length;
		
		chip = new ChiSquaredDistributionImpl(test*2);

		int cnt = 0;
		for (int i = 0; i < phe1.length; i++)
		{
			for (int j = 0; j < phe2.length; j++)
			{
				double[] s = new double[test];
				double[] p = new double[test];
				for (int k = 0; k < test; k++)
				{
					s[k] = Math.abs(phe1[i][k] - phe2[j][k]);
					try
					{
						p[k] = (ndp.cumulativeProbability(s[k]) - 0.5)*2;
						if (p[k] == 0)
						{
							p[k] = 10e-20;
						}
					}
					catch (MathException e)
					{
						e.printStackTrace();
					}
				}

				cnt++;
				predictorFile.print(ID1.get(i) + "\t" + ID2.get(j) + "\t");
				double logp = 0;
				double chi = 0;
				for (int l = 0; l < test; l++)
				{
					logp += -1 * Math.log10(p[l]);
					chi += -1 * Math.log(p[l]);
					if (l < (test - 1))
					{
						predictorFile.print(s[l] + "\t" + p[l] + "\t");
					}
					else
					{
						double cp = 0;
						try
						{
							cp = 1 - chip.cumulativeProbability(chi);
						}
						catch (MathException e)
						{
							e.printStackTrace();
						}
						predictorFile.println(s[l] + "\t" + p[l] + "\t" + logp + "\t" + chi + "\t" + cp);
					}
				}
			}
		}
		predictorFile.close();
		Logger.printUserLog("In total " + cnt + " similar pairs were detected.");
	}

	private void readPhenotypes1(String fileName)
	{
		InputDataSet ds = new InputDataSet();
		Logger.printUserLog("Reading " + fileName);
		ds.readPhenotypeFile(fileName);
//		gear.util.BufferedReader reader = gear.util.BufferedReader
//				.openTextFile(fileName, "phenotype");
//		HashSet<SubjectID> subjectsRead = new HashSet<SubjectID>();
//		String[] tokens = reader.readTokensAtLeast(3);
//		if (tokens == null)
//		{
//			Logger.printUserError("The phenotype file '" + fileName + "' is empty.");
//			System.exit(1);
//		}
//		int numCols = tokens.length;
//		ArrayList<String> tmp = NewIt.newArrayList();
//		int indIdx = 0;
//		while ((tokens = reader.readTokens(numCols)) != null)
//		{
//			SubjectID subID = new SubjectID(/* famID */tokens[0], /* indID */
//					tokens[1]);
//			if (subjectsRead.contains(subID))
//			{
//				Logger.printUserLog("Individual:" + subID.getFamilyID() + " " + subID
//						.getIndividualID() + " at line " + (indIdx + 1) + " excluded due to double entry.");
//			}
//			else
//			{
//				StringBuffer id = new StringBuffer();
//				id.append(subID.getFamilyID() + "\t" + subID.getIndividualID());
//				StringBuffer sb = new StringBuffer();
//				sb.append(id.toString());
//				for (int i = 2; i < numCols; i++)
//				{
//					sb.append("\t" + tokens[i]);
//				}
//				tmp.add(sb.toString());
//			}
//			indIdx++;
//		}
//		reader.close();
		phe1 = new double[ds.getNumberOfSubjects()][ds.getNumberOfTraits()];
		for (int i = 0; i < ds.getNumberOfSubjects(); i++)
		{
			StringBuffer id = new StringBuffer();
			ID1.add(ds.getSubjectID(i).toString());
			for (int j = 0; j < ds.getNumberOfTraits(); j++)
			{
				phe1[i][j] = ds.getPhenotype(i, j);
			}
		}
		Logger.printUserLog("Read " + ds.getNumberOfSubjects() + " individuals, " + ds.getNumberOfTraits() + " scores. in " + fileName);
//		phe1 = NCGStatUtils.standardize(phe1, true);
	}

	private void readPhenotypes2(String fileName)
	{
		InputDataSet ds = new InputDataSet();
		Logger.printUserLog("Reading " + fileName);
		ds.readPhenotypeFile(fileName);

//		Logger.printUserLog("Reading " + fileName);
//		gear.util.BufferedReader reader = gear.util.BufferedReader
//				.openTextFile(fileName, "phenotype");
//		HashSet<SubjectID> subjectsRead = new HashSet<SubjectID>();
//		String[] tokens = reader.readTokensAtLeast(3);
//		if (tokens == null)
//		{
//			Logger.printUserError("The phenotype file '" + fileName + "' is empty.");
//			System.exit(1);
//		}
//		int numCols = tokens.length;
//		ArrayList<String> tmp = NewIt.newArrayList();
//		int indIdx = 0;
//		while ((tokens = reader.readTokens(numCols)) != null)
//		{
//			SubjectID subID = new SubjectID(/* famID */tokens[0], /* indID */
//					tokens[1]);
//			if (subjectsRead.contains(subID))
//			{
//				Logger.printUserLog("Individual:" + subID.getFamilyID() + " " + subID
//						.getIndividualID() + " at line " + (indIdx + 1) + " excluded due to double entry.");
//			}
//			else
//			{
//				StringBuffer id = new StringBuffer();
//				id.append(subID.getFamilyID() + "\t" + subID.getIndividualID());
//				StringBuffer sb = new StringBuffer();
//				sb.append(id.toString());
//				for (int i = 2; i < numCols; i++)
//				{
//					sb.append("\t" + tokens[i]);
//				}
//				tmp.add(sb.toString());
//			}
//			indIdx++;
//		}
//		reader.close();
		phe2 = new double[ds.getNumberOfSubjects()][ds.getNumberOfTraits()];
		for (int i = 0; i < ds.getNumberOfSubjects(); i++)
		{
			ID2.add(ds.getSubjectID(i).toString());
			for (int j = 0; j < ds.getNumberOfTraits(); j++)
			{
				phe2[i][j] = ds.getPhenotype(i, j);
			}
		}
		Logger.printUserLog("Read " + ds.getNumberOfSubjects() + " individuals, " + ds.getNumberOfTraits() + " scores. in " + fileName);
//		phe2 = NCGStatUtils.standardize(phe2, true);
	}
}
