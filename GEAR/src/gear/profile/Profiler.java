package gear.profile;

import gear.CmdArgs;
import gear.ConstValues;
import gear.family.pedigree.file.SNP;
import gear.util.BufferedReader;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.SNPMatch;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;

public class Profiler
{	
	public static void makeProfile()
	{
		checkCmdArgs();
		
		HashMap<String, Score> scoreMap = readScores();  // LocusName-to-Score map
		HashMap<String, Float> qScoreMap = readQScores();  // LocusName-to-QScore map
		QRange[] qRanges = readQRanges();
		
		int[] matchNums = new int[AlleleMatchScheme.values().length];
		
		Data genoData = initData();
		SNP[] snps = genoData.getSNPs();
		
		int numLociNoScore = 0, numMonoLoci = 0, numAmbiguousLoci = 0, numLociNoQScore = 0;
		Score[] scores = new Score[snps.length]; 
		
		AlleleMatchScheme[] alleleMatchSchemes = new AlleleMatchScheme[snps.length];
		for (int i = 0; i < alleleMatchSchemes.length; ++i)
		{
			alleleMatchSchemes[i] = AlleleMatchScheme.MATCH_NONE;
		}
		
		int numLocusGroups = qRanges == null ? 1 : qRanges.length;
		
		boolean[][] isInLocusGroup = new boolean[snps.length][numLocusGroups];
		int[] numInLocusGroup = new int[numLocusGroups];
		
		// Check whether each SNP should be used for profiling.
		for (int snpIdx = 0; snpIdx < snps.length; ++snpIdx)
		{
			SNP snp = snps[snpIdx];
			
			Score score = scoreMap.remove(snp.getName());
			scores[snpIdx] = score;
			
			if (score == null)
			{
				++numLociNoScore;
				continue;
			}
			
			if (CmdArgs.INSTANCE.getTranFunction() == gear.RegressionModel.LOGIT)
			{
				scores[snpIdx].setValue((float)Math.log(score.getValue()));
			}
			
			char allele1 = snp.getFirstAllele();
			char allele2 = snp.getSecAllele();
			
			if (!SNPMatch.isBiallelic(allele1, allele2))
			{
				++numMonoLoci;
				continue;
			}
			
			if (SNPMatch.isAmbiguous(allele1, allele2))  // A/T or C/G
			{
				++numAmbiguousLoci;
				if (!CmdArgs.INSTANCE.keepATGC())
				{
					continue;
				}
			}
			
			AlleleMatchScheme matchScheme = getMatchScheme(score.getAllele(), allele1, allele2);
			++matchNums[matchScheme.ordinal()];
			if (matchScheme == AlleleMatchScheme.MATCH_NONE)
			{
				continue;
			}
			alleleMatchSchemes[snpIdx] = matchScheme;
			
			if (qScoreMap == null)
			{
				isInLocusGroup[snpIdx][0] = true;
			}
			else
			{
				Float qScore = qScoreMap.remove(snp.getName());
				
				if (qScore == null)
				{
					++numLociNoQScore;
					continue;
				}

				assert qRanges != null;

				for (int rangeIdx = 0; rangeIdx < qRanges.length; ++rangeIdx)
				{
					QRange qRange = qRanges[rangeIdx];
					isInLocusGroup[snpIdx][rangeIdx] = qRange.getLowerBound() <= qScore && qScore <= qRange.getUpperBound();
					if (isInLocusGroup[snpIdx][rangeIdx])
					{
						numInLocusGroup[rangeIdx]++;
					}
				}
			}
		}  // for each SNP

		Logger.printUserLog("Number of loci having no score (because they do not appear in the score file, or their scores are invalid, etc.): " + numLociNoScore);
		Logger.printUserLog("Number of monomorphic loci (removed): " + numMonoLoci);
		Logger.printUserLog("Number of ambiguous loci (A/T or C/G) " + (CmdArgs.INSTANCE.keepATGC() ? "detected: " : "removed: ") + numAmbiguousLoci);

		// Allele Matching Schemes
		Logger.printUserLog("Number of Scheme I predictors: predictor alleles were A1: " + matchNums[AlleleMatchScheme.MATCH_ALLELE1.ordinal()]);
		Logger.printUserLog("Number of Scheme II predictors: predictor alleles were A2: " + matchNums[AlleleMatchScheme.MATCH_ALLELE2.ordinal()]);
		if (!CmdArgs.INSTANCE.getProfileArgs().isAutoFlipOff())
		{
			Logger.printUserLog("Number of Scheme III predictors: predictor alleles were flipped A1: " + matchNums[AlleleMatchScheme.MATCH_ALLELE1_FLIPPED.ordinal()]);
			Logger.printUserLog("Number of Scheme IV predictors: predictor alleles were flipped A2: " + matchNums[AlleleMatchScheme.MATCH_ALLELE2_FLIPPED.ordinal()]);
		}
		Logger.printUserLog("Number of score alleles matching none in the data file: " + matchNums[AlleleMatchScheme.MATCH_NONE.ordinal()]);

		if (qScoreMap != null)
		{
			Logger.printUserLog("Number of loci having no q-scores: " + numLociNoQScore);
			for (int i = 0; i < numInLocusGroup.length; i++)
			{
				QRange qRange = qRanges[i];
				Logger.printUserLog("\tNumber of loci within the range: " + qRange.getLowerBound() + ", " + qRange.getUpperBound() + " is " + numInLocusGroup[i]);
			}
		}

		ArrayList<String> famIDs = new ArrayList<String>();
		ArrayList<String> indIDs = new ArrayList<String>();
		ArrayList<String> phenos = new ArrayList<String>();
		ArrayList<float[]> riskProfiles = new ArrayList<float[]>();
		ArrayList<int[]> numLociUsed = new ArrayList<int[]>();
		
		Data.Iterator iter = genoData.iterator();
		while (iter.next())
		{
			int indIdx = iter.getIndividualIndex();
			int locIdx = iter.getLocusIndex();
			
			while (indIdx >= famIDs.size())
			{
				famIDs.add(null);
			}
			famIDs.set(indIdx, iter.getFamilyID());
			
			while (indIdx >= indIDs.size())
			{
				indIDs.add(null);
			}
			indIDs.set(indIdx, iter.getIndividualID());
			
			while (indIdx >= phenos.size())
			{
				phenos.add(null);
			}
			phenos.set(indIdx, iter.getPhenotype());
			
			while (indIdx >= riskProfiles.size())
			{
				riskProfiles.add(new float [numLocusGroups]);
			}
			
			while (indIdx >= numLociUsed.size())
			{
				numLociUsed.add(new int [numLocusGroups]);
			}
			
			if (scores[locIdx] == null)
			{
				continue;
			}
			
			float scoreAlleleFrac = 0.0f;
			switch (alleleMatchSchemes[locIdx])
			{
			case MATCH_ALLELE1:
			case MATCH_ALLELE1_FLIPPED:
				scoreAlleleFrac = iter.getAllele1Fraction();
				break;
			case MATCH_ALLELE2:
			case MATCH_ALLELE2_FLIPPED:
				scoreAlleleFrac = 2.0f - iter.getAllele1Fraction();
				break;
			default:
				continue;
			}
			
			float riskValue = coeffModel.compute(scoreAlleleFrac) * scores[locIdx].getValue();
			
			for (int locGrpIdx = 0; locGrpIdx < isInLocusGroup[locIdx].length; ++locGrpIdx)
			{
				if (isInLocusGroup[locIdx][locGrpIdx])
				{
					riskProfiles.get(indIdx)[locGrpIdx] += riskValue;
					++numLociUsed.get(indIdx)[locGrpIdx];
				}
			}
		}

		for (int indIdx = 0; indIdx < riskProfiles.size(); ++indIdx)
		{
			for (int locGrpIdx = 0; locGrpIdx < numLocusGroups; ++locGrpIdx)
			{
				int denom = numLociUsed.get(indIdx)[locGrpIdx];
				if (denom != 0)
				{
					riskProfiles.get(indIdx)[locGrpIdx] /= sameAsPlink ? (denom << 1) : denom;
				}
			}
		}
		
		PrintStream predictorFile = FileUtil.CreatePrintStream(resultFile);
		
		// Title Line
		predictorFile.print("FID\tIID\tPHENO");
		if (qRanges == null)
		{
			predictorFile.print("\tSCORE");
		}
		else
		{
			for (int rangeIdx = 0; rangeIdx < qRanges.length; ++rangeIdx)
			{
				predictorFile.print("\tSCORE." + qRanges[rangeIdx].getName());
			}
		}
		predictorFile.println();
		
		for (int indIdx = 0; indIdx < riskProfiles.size(); indIdx++)
		{
			predictorFile.print(famIDs.get(indIdx) + "\t" + indIDs.get(indIdx) + "\t" + phenos.get(indIdx));
			for (int locGrpIdx = 0; locGrpIdx < numLocusGroups; ++locGrpIdx)
			{
				predictorFile.print("\t" + riskProfiles.get(indIdx)[locGrpIdx]);
			}
			predictorFile.println();
		}
		predictorFile.close();
	}

	private static void checkCmdArgs()
	{
		String scoreFile = CmdArgs.INSTANCE.getProfileArgs().getScoreFile();
		if (scoreFile == null)
		{
			Logger.printUserError("Score file is not provided.");
			System.exit(1);
		}
		
		String qScoreFile = CmdArgs.INSTANCE.getProfileArgs().getQScoreFile();
		String qRangeFile = CmdArgs.INSTANCE.getProfileArgs().getQRangeFile();
		if (qScoreFile == null && qRangeFile != null || 
			qScoreFile != null && qRangeFile == null)
		{
			Logger.printUserError("--q-score-file and --q-score-range must be set together.");
			System.exit(1);
		}
		
		if (!CmdArgs.INSTANCE.getFileArgs().isSet() && !CmdArgs.INSTANCE.getBFileArgs(0).isSet())
		{
			String dosage = CmdArgs.INSTANCE.getProfileArgs().getMachDosageFile();
			String info = CmdArgs.INSTANCE.getProfileArgs().getMachInfoFile();
			String dosageBatch = CmdArgs.INSTANCE.getProfileArgs().getMachDosageBatchFile();
			String infoBatch = CmdArgs.INSTANCE.getProfileArgs().getMachInfoBatchFile();
			
			if (dosage == null && info == null && dosageBatch == null && infoBatch == null)
			{
				Logger.printUserError("No input data files.");
				System.exit(1);
			}
			
			if (dosage == null && info != null || dosage != null && info == null)
			{
				Logger.printUserError("--mach-dosage and --mach-info must be set together.");
				System.exit(1);
			}
			
			if (dosageBatch == null && infoBatch != null || dosageBatch != null && infoBatch == null)
			{
				Logger.printUserError("--mach-dosage-batch and --mach-info-batch must be set together.");
				System.exit(1);
			}
		}
		
		resultFile = CmdArgs.INSTANCE.out;
		
		// Compute Model
		String sCoeffModel = CmdArgs.INSTANCE.getProfileArgs().getModel();
		if (sCoeffModel == null)
		{
			Logger.printUserLog("Allelic model is used.");
			coeffModel = new AdditiveCoeffModel();
			sameAsPlink = true;
		}
		else
		{
			resultFile += "." + sCoeffModel;
			if (sCoeffModel.equals("additive"))
			{
				Logger.printUserLog("Additive model is used.");
				coeffModel = new AdditiveCoeffModel();
			}
			else if (sCoeffModel.equals("dominance"))
			{
				Logger.printUserLog("Dominance model is used.");
				coeffModel = new DominanceCoeffModel();
			}
			else if (sCoeffModel.equals("recessive"))
			{
				Logger.printUserLog("Recessive model is used.");
				coeffModel = new RecessiveCoeffModel();
			}
			else
			{
				String msg = "";
				msg += "'" + sCoeffModel + "' is an invalid coefficient model. ";
				msg += "Valid models are 'additive', 'dominance' or 'recessive'.";
				Logger.printUserError(msg);
				System.exit(1);
			}
		}
		
		resultFile += ".profile";
	}
	
	private static HashMap<String, Score> readScores()
	{
		HashMap<String, Score> scores = new HashMap<String, Score>();
		
		BufferedReader reader = BufferedReader.openTextFile(CmdArgs.INSTANCE.getProfileArgs().getScoreFile(), "score");
		String[] tokens;
		
		while ((tokens = reader.readTokens(3)) != null)
		{
			if (tokens[1].length() != 1)
			{
				reader.errorPreviousLine("The allele is not a character.");
			}
			
			if (!ConstValues.isNA(tokens[2]))
			{
				try
				{
					scores.put(/* locusName = */ tokens[0], new Score(/* scoreAllele = */ tokens[1].charAt(0), /* score = */ Float.parseFloat(tokens[2])));
				}
				catch (NumberFormatException e)
				{
					reader.errorPreviousLine("'" + tokens[2] + "' is not a valid score.");
				}
			}
		}
		reader.close();
		
		Logger.printUserLog("Number of valid scores: " + scores.size());
		
		return scores;
	}

	private static HashMap<String, Float> readQScores()
	{
		String qScoreFile = CmdArgs.INSTANCE.getProfileArgs().getQScoreFile();
		if (qScoreFile == null)
		{
			return null;
		}
		
		HashMap<String, Float> qScores = new HashMap<String, Float>();

		BufferedReader reader = BufferedReader.openTextFile(qScoreFile, "q-score");
		String[] tokens;
		
		while ((tokens = reader.readTokens(2)) != null)
		{
			if (!ConstValues.isNA(tokens[1]))
			{
				try
				{
					qScores.put(/* locusName = */ tokens[0], Float.parseFloat(tokens[1]));
				}
				catch (NumberFormatException e)
				{
					reader.errorPreviousLine("'" + tokens[1] + "' is not a valid q-score.");
				}
			}
		}
		reader.close();

		Logger.printUserLog("Number of q-scores: " + qScores.size());
		
		return qScores;
	}

	private static QRange[] readQRanges()
	{
		String qRangeFile = CmdArgs.INSTANCE.getProfileArgs().getQRangeFile();
		if (qRangeFile == null)
		{
			return null;
		}
		
		ArrayList<QRange> qRanges = new ArrayList<QRange>();

		BufferedReader reader = BufferedReader.openTextFile(qRangeFile, "q-score-range");
		String[] tokens;
		
		while ((tokens = reader.readTokens(3)) != null)
		{
			try
			{
				float lowerBound, upperBound;
				lowerBound = ConstValues.isNA(tokens[1]) ? Float.MIN_VALUE : Float.parseFloat(tokens[1]);
				upperBound = ConstValues.isNA(tokens[2]) ? Float.MAX_VALUE : Float.parseFloat(tokens[2]);
				qRanges.add(new QRange(/* name = */tokens[0], lowerBound, upperBound));
			}
			catch (NumberFormatException e)
			{
				reader.errorPreviousLine("An invalid q-score value is detected.");
			}
		}
		reader.close();

		Logger.printUserLog("Number of q-ranges: " + qRanges.size());

		return qRanges.toArray(new QRange[0]);
	}
	
	private static Data initData()
	{
		if (CmdArgs.INSTANCE.getFileArgs().isSet() || CmdArgs.INSTANCE.getBFileArgs(0).isSet())
		{
			return PlinkData.create();
		}
		return MachData.create();
	}
	
	private static AlleleMatchScheme getMatchScheme(char scoreAllele, char allele1, char allele2)
	{
		if (scoreAllele == allele1)
		{
			return AlleleMatchScheme.MATCH_ALLELE1;
		}
		else if (scoreAllele == allele2)
		{
			return AlleleMatchScheme.MATCH_ALLELE2;
		}
		else if (!CmdArgs.INSTANCE.getProfileArgs().isAutoFlipOff())
		{
			if (scoreAllele == SNPMatch.Flip(allele1))
			{
				return AlleleMatchScheme.MATCH_ALLELE1_FLIPPED;
			}
			else if (scoreAllele == SNPMatch.Flip(allele2))
			{
				return AlleleMatchScheme.MATCH_ALLELE2_FLIPPED;
			}
		}
		return AlleleMatchScheme.MATCH_NONE;
	}

	private static CoeffModel coeffModel;
	private static boolean sameAsPlink;
	private static String resultFile;
}
