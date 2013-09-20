package gear.profile;

import gear.CommandArguments;
import gear.CommandImpl;
import gear.ConstValues;
import gear.family.pedigree.file.SNP;
import gear.util.BufferedReader;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.SNPMatch;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;

public final class ProfileCommandImpl extends CommandImpl
{
	@Override
	public void execute(CommandArguments cmdArgs)
	{
		profCmdArgs = (ProfileCommandArguments)cmdArgs;
		
		printWhichCoeffModelIsUsed();
		
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
			
			if (profCmdArgs.getIsLogit())
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
				if (!profCmdArgs.getIsKeepATGC())
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
		Logger.printUserLog("Number of ambiguous loci (A/T or C/G) " + (profCmdArgs.getIsKeepATGC() ? "detected: " : "removed: ") + numAmbiguousLoci);

		// Allele Matching Schemes
		Logger.printUserLog("Number of Scheme I predictors: predictor alleles were A1: " + matchNums[AlleleMatchScheme.MATCH_ALLELE1.ordinal()]);
		Logger.printUserLog("Number of Scheme II predictors: predictor alleles were A2: " + matchNums[AlleleMatchScheme.MATCH_ALLELE2.ordinal()]);
		if (profCmdArgs.getIsAutoFlip())
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
			
			float riskValue = profCmdArgs.getCoeffModel().compute(scoreAlleleFrac) * scores[locIdx].getValue();
			
			for (int locGrpIdx = 0; locGrpIdx < isInLocusGroup[locIdx].length; ++locGrpIdx)
			{
				if (isInLocusGroup[locIdx][locGrpIdx])
				{
					riskProfiles.get(indIdx)[locGrpIdx] += riskValue;
					++numLociUsed.get(indIdx)[locGrpIdx];
				}
			}
		}

		if (profCmdArgs.getIsWeighted())
		{
			for (int indIdx = 0; indIdx < riskProfiles.size(); ++indIdx)
			{
				for (int locGrpIdx = 0; locGrpIdx < numLocusGroups; ++locGrpIdx)
				{
					int denom = numLociUsed.get(indIdx)[locGrpIdx];
					if (denom != 0)
					{
						riskProfiles.get(indIdx)[locGrpIdx] /= profCmdArgs.getIsSameAsPlink() ? (denom << 1) : denom;
					}
				}
			}
		}
		PrintStream predictorFile = FileUtil.CreatePrintStream(profCmdArgs.getResultFile());
		
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
	
	private void printWhichCoeffModelIsUsed()
	{
		if (profCmdArgs.getIsSameAsPlink())
		{
			Logger.printUserLog("PLINK allelic model is used.");
		}
		else if (profCmdArgs.getCoeffModel() instanceof AdditiveCoeffModel)
		{
			Logger.printUserLog("Additive model is used.");
		}
		else if (profCmdArgs.getCoeffModel() instanceof DominanceCoeffModel)
		{
			Logger.printUserLog("Dominance model is used.");
		}
		else
		{
			Logger.printUserLog("Recessive model is used.");
		}
	}
	
	private HashMap<String, Score> readScores()
	{
		HashMap<String, Score> scores = new HashMap<String, Score>();
		
		BufferedReader reader = BufferedReader.openTextFile(profCmdArgs.getScoreFile(), "score");
		String[] tokens;

		while ((tokens = reader.readTokens(3)) != null)
		{
			if (tokens[1].length() != 1)
			{
				reader.warningPreviousLine("'" + tokens[1] + "' is not a character, so it is not a valid allele, and this line will be ignored.");
				continue;
			}
			
			if (!ConstValues.isNA(tokens[2]))
			{
				try
				{
					scores.put(/* locusName = */ tokens[0], new Score(/* scoreAllele = */ tokens[1].charAt(0), /* score = */ Float.parseFloat(tokens[2])));
				}
				catch (NumberFormatException e)
				{
					reader.warningPreviousLine("'" + tokens[2] + "' is not a floating point number, so it it not a valid score, and this line will be ingored.");
				}
			}
		}
		reader.close();
		
		Logger.printUserLog("Number of valid scores: " + scores.size());
		return scores;
	}

	private HashMap<String, Float> readQScores()
	{
		String qScoreFile = profCmdArgs.getQScoreFile();
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

	private QRange[] readQRanges()
	{
		String qRangeFile = profCmdArgs.getQRangeFile();
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
	
	private Data initData()
	{
		if (profCmdArgs.getFile() != null)
		{
			return PlinkData.createByFile(profCmdArgs.getFile());
		}
		else if (profCmdArgs.getBFile() != null)
		{
			return PlinkData.createByBFile(profCmdArgs.getBFile());
		}
		return MachData.create(profCmdArgs.getMachDosageFile(),
		                       profCmdArgs.getMachInfoFile(),
		                       profCmdArgs.getMachDosageBatch(),
		                       profCmdArgs.getMachInfoBatch());
	}
	
	private AlleleMatchScheme getMatchScheme(char scoreAllele, char allele1, char allele2)
	{
		if (scoreAllele == allele1)
		{
			return AlleleMatchScheme.MATCH_ALLELE1;
		}
		else if (scoreAllele == allele2)
		{
			return AlleleMatchScheme.MATCH_ALLELE2;
		}
		else if (profCmdArgs.getIsAutoFlip())
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

	private ProfileCommandArguments profCmdArgs;
}
