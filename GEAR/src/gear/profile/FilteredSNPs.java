package gear.profile;

import gear.family.pedigree.file.SNP;
import gear.util.SNPMatch;

import java.util.HashMap;

class FilteredSNPs
{
	private FilteredSNPs(int numSNPs, int numLocusGroups)
	{
		scores = new Score[numSNPs];
		matchNums = new int[AlleleMatchScheme.values().length];
		alleleMatchSchemes = new AlleleMatchScheme[numSNPs];
		for (int i = 0; i < alleleMatchSchemes.length; ++i)
		{
			alleleMatchSchemes[i] = AlleleMatchScheme.MATCH_NONE;
		}
		isInLocusGroup = new boolean[numSNPs][numLocusGroups];
		numInLocusGroup = new int[numLocusGroups];
	}
	
	// Check whether each SNP should be used for profiling.
	public static FilteredSNPs filter(	SNP[]					snps,
										HashMap<String, Score>	scoreMap,
										HashMap<String, Float>	qScoreMap,
										QRange[]				qRanges,
										ProfileCommandArguments	cmdArgs	)
	{
		FilteredSNPs ret = new FilteredSNPs(snps.length, /*numLocusGroups*/qRanges == null ? 1 : qRanges.length);
		
		for (int snpIdx = 0; snpIdx < snps.length; ++snpIdx)
		{
			SNP snp = snps[snpIdx];
			
			Score score = scoreMap.remove(snp.getName());
			ret.scores[snpIdx] = score;
			
			if (score == null)
			{
				++ret.numLociNoScore;
				continue;
			}
			
			if (cmdArgs.getIsLogit())
			{
				ret.scores[snpIdx].setValue((float)Math.log(score.getValue()));
			}
			
			char allele1 = snp.getFirstAllele();
			char allele2 = snp.getSecAllele();
			
			if (!SNPMatch.isBiallelic(allele1, allele2))
			{
				++ret.numMonoLoci;
				continue;
			}
			
			if (SNPMatch.isAmbiguous(allele1, allele2))  // A/T or C/G
			{
				++ret.numAmbiguousLoci;
				if (!cmdArgs.getIsKeepATGC())
				{
					continue;
				}
			}
			
			AlleleMatchScheme matchScheme = getMatchScheme(score.getAllele(), allele1, allele2, cmdArgs.getIsAutoFlip());
			++ret.matchNums[matchScheme.ordinal()];
			if (matchScheme == AlleleMatchScheme.MATCH_NONE)
			{
				continue;
			}
			ret.alleleMatchSchemes[snpIdx] = matchScheme;
			
			if (qScoreMap == null)
			{
				ret.isInLocusGroup[snpIdx][0] = true;
			}
			else
			{
				Float qScore = qScoreMap.remove(snp.getName());
				
				if (qScore == null)
				{
					++ret.numLociNoQScore;
					continue;
				}

				assert qRanges != null;

				for (int rangeIdx = 0; rangeIdx < qRanges.length; ++rangeIdx)
				{
					QRange qRange = qRanges[rangeIdx];
					ret.isInLocusGroup[snpIdx][rangeIdx] = qRange.getLowerBound() <= qScore && qScore <= qRange.getUpperBound();
					if (ret.isInLocusGroup[snpIdx][rangeIdx])
					{
						ret.numInLocusGroup[rangeIdx]++;
					}
				}
			}
		}  // for each SNP
		
		return ret;
	}
	
	private static AlleleMatchScheme getMatchScheme(char scoreAllele, char allele1, char allele2, boolean isAutoFlip)
	{
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
	
	public Score getScore(int locIdx)
	{
		return scores[locIdx];
	}
	
	public int getNumLociNoScore()
	{
		return numLociNoScore;
	}
	
	public int getNumMonoLoci()
	{
		return numMonoLoci;
	}
	
	public int getNumAmbiguousLoci()
	{
		return numAmbiguousLoci;
	}
	
	public int getNumLociNoQScore()
	{
		return numLociNoQScore;
	}
	
	public int getMatchNum(AlleleMatchScheme scheme)
	{
		return matchNums[scheme.ordinal()];
	}
	
	public AlleleMatchScheme getAlleleMatchScheme(int locIdx)
	{
		return alleleMatchSchemes[locIdx];
	}
	
	public boolean isInLocusGroup(int locIdx, int locGrpIdx)
	{
		return isInLocusGroup[locIdx][locGrpIdx];
	}
	
	public int getNumLocusGroups()
	{
		return numInLocusGroup.length;
	}
	
	public int getNumInLocusGroup(int locGrpIdx)
	{
		return numInLocusGroup[locGrpIdx];
	}
	
	private Score[] scores;
	private int numLociNoScore = 0, numMonoLoci = 0, numAmbiguousLoci = 0, numLociNoQScore = 0;
	private int[] matchNums;
	private AlleleMatchScheme[] alleleMatchSchemes;
	private boolean[][] isInLocusGroup;
	private int[] numInLocusGroup;
}
