package gear.subcommands.profile;

import gear.family.pedigree.file.SNP;
import gear.util.SNPMatch;

import java.util.HashMap;

class FilteredSNPs
{
	protected FilteredSNPs(	SNP[]					snps,
							ScoreFile				scoreFile,
							HashMap<String, Float>	qScoreMap,
							QRange[]				qRanges,
							ProfileCommandArguments	cmdArgs	)
	{
		init(scoreFile.getNumberOfTraits(), snps.length, /*numLocusGroups*/qRanges == null ? 1 : qRanges.length);
		
		Score[] scores = arrangeScores(snps, scoreFile);
		
		convertScores(scores, scoreFile.getNumberOfTraits(), cmdArgs.getIsLogit());

		for (int snpIdx = 0; snpIdx < snps.length; ++snpIdx)
		{
			SNP snp = snps[snpIdx];
			
			if (scores[snpIdx] == null)
			{
				continue;
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
				if (!cmdArgs.getIsKeepATGC())
				{
					continue;
				}
			}
			
			AlleleMatchScheme matchScheme = getMatchScheme(scores[snpIdx].getAllele(), allele1, allele2, cmdArgs.getIsAutoFlip());
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
				Float qScore = qScoreMap.get(snp.getName());
				
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
	}
	
	private void init(int numTraits, int numSNPs, int numLocusGroups)
	{
		scores = new Float[numTraits][numSNPs];
		matchNums = new int[AlleleMatchScheme.values().length];
		alleleMatchSchemes = new AlleleMatchScheme[numSNPs];
		for (int i = 0; i < alleleMatchSchemes.length; ++i)
		{
			alleleMatchSchemes[i] = AlleleMatchScheme.MATCH_NONE;
		}
		isInLocusGroup = new boolean[numSNPs][numLocusGroups];
		numInLocusGroup = new int[numLocusGroups];
	}
	
	private static Score[] arrangeScores(SNP[] snps, ScoreFile scoreFile)
	{
		Score[] scores = new Score[snps.length];	
		for (int snpIdx = 0; snpIdx < snps.length; ++snpIdx)
		{
			scores[snpIdx] = scoreFile.getScore(snps[snpIdx].getName());
		}
		return scores;
	}
	
	private void convertScores(Score[] scores, int numTraits, boolean isLogit)
	{
		this.scores = new Float[numTraits][scores.length];
		
		for (int traitIdx = 0; traitIdx < numTraits; ++traitIdx)
		{
			for (int snpIdx = 0; snpIdx < scores.length; ++snpIdx)
			{
				Score score = scores[snpIdx]; 
				if (score != null && score.getValue(traitIdx) != null)
				{
					if (isLogit)
					{
						if (score.getValue(traitIdx).floatValue() > 0)
						{
							this.scores[traitIdx][snpIdx] = (float)Math.log(score.getValue(traitIdx));
						}
					}
					else
					{
						this.scores[traitIdx][snpIdx] = (float) score.getValue(traitIdx);
					}
				}
			}
		}
	}
	
	protected static AlleleMatchScheme getMatchScheme(char scoreAllele, char allele1, char allele2, boolean isAutoFlip)
	{
		if (scoreAllele >= 97 && scoreAllele <= 122)
		{
			scoreAllele -= 32;
		}
		
		if (allele1 >= 97 && allele1 <= 122)
		{
			allele1 -= 32;
		}

		if (allele2 >= 97 && allele2 <= 122)
		{
			allele2 -= 32;
		}

		if (scoreAllele == allele1)
		{
			return AlleleMatchScheme.MATCH_ALLELE1;
		}
		else if (scoreAllele == allele2)
		{
			return AlleleMatchScheme.MATCH_ALLELE2;
		}
		else if (isAutoFlip)
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
	
	protected int getNumLociNoScore(int traitIdx)
	{
		int ret = 0;
		for (Float score : scores[traitIdx])
		{
			if (score == null)
			{
				++ret;
			}
		}
		return ret;
	}
	
	protected Float getScore(int traitIdx, int snpIdx)
	{
		return scores[traitIdx][snpIdx];
	}
	
	protected int getNumMonoLoci()
	{
		return numMonoLoci;
	}
	
	protected int getNumAmbiguousLoci()
	{
		return numAmbiguousLoci;
	}
	
	protected int getNumLociNoQScore()
	{
		return numLociNoQScore;
	}
	
	protected int getMatchNum(AlleleMatchScheme scheme)
	{
		return matchNums[scheme.ordinal()];
	}
	
	protected AlleleMatchScheme getAlleleMatchScheme(int locIdx)
	{
		return alleleMatchSchemes[locIdx];
	}
	
	protected boolean isInLocusGroup(int locIdx, int locGrpIdx)
	{
		return isInLocusGroup[locIdx][locGrpIdx];
	}
	
	protected int getNumLocusGroups()
	{
		return numInLocusGroup.length;
	}
	
	protected int getNumInLocusGroup(int locGrpIdx)
	{
		return numInLocusGroup[locGrpIdx];
	}
	
	private Float[][] scores;
	private int numMonoLoci = 0, numAmbiguousLoci = 0, numLociNoQScore = 0;
	private int[] matchNums;
	private AlleleMatchScheme[] alleleMatchSchemes;
	private boolean[][] isInLocusGroup;
	private int[] numInLocusGroup;
}
