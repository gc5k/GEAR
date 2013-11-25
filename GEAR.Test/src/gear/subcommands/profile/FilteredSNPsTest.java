package gear.subcommands.profile;

import static org.junit.Assert.*;

import java.util.HashMap;

import gear.ConstValues;
import gear.family.pedigree.file.SNP;
import gear.subcommands.profile.AlleleMatchScheme;
import gear.subcommands.profile.FilteredSNPs;
import gear.subcommands.profile.ProfileCommandArguments;
import gear.subcommands.profile.QRange;
import gear.subcommands.profile.ScoreFile;

import org.junit.Test;

public class FilteredSNPsTest
{

	@Test
	public void test()
	{
		SNP[] snps = new SNP[10];
		snps[0] = new SNP("SNP1", 'A', 'C');
		snps[1] = new SNP("SNP2", 'A', 'T');
		snps[2] = new SNP("SNP3", 'A', 'G');
		snps[3] = new SNP("SNP4", 'A', ConstValues.MISSING_ALLELE_CHAR);
		snps[4] = new SNP("SNP5", 'C', 'A');
		snps[5] = new SNP("SNP6", 'C', 'T');
		snps[6] = new SNP("SNP7", 'C', 'G');
		snps[7] = new SNP("SNP8", ConstValues.MISSING_ALLELE_CHAR, 'C');
		snps[8] = new SNP("SNP9", 'T', 'C');
		snps[9] = new SNP("SNP10", 'T', 'G');
		
		ScoreFile scoreFile = new ScoreFile("data/ScoresWithHeaders.txt", /*hasHeaders*/true);
		
		HashMap<String, Float> qScoreMap = new HashMap<String, Float>();
		qScoreMap.put("SNP1", 2.0f);
		qScoreMap.put("SNP2", 3.0f);
		qScoreMap.put("SNP3", 3.5f);
		qScoreMap.put("SNP4", 2.2f);
		qScoreMap.put("SNP5", 4.1f);
		qScoreMap.put("SNP6", 4.4f);
		qScoreMap.put("SNP10", 3.3f);
		
		QRange[] qRanges = new QRange[4];
		qRanges[0] = new QRange("QRange1", 1.5f, 2.5f);
		qRanges[1] = new QRange("QRange2", 2.1f, 3.4f);
		qRanges[2] = new QRange("QRange3", 3.2f, 4.0f);
		qRanges[3] = new QRange("QRange4", 4.0f, 5.0f);
		
		ProfileCommandArguments cmdArgs = new ProfileCommandArguments();
		cmdArgs.setIsKeepATGC(true);
		cmdArgs.setIsAutoFlip(true);
		
		FilteredSNPs filteredSNPs = new FilteredSNPs(snps, scoreFile, qScoreMap, qRanges, cmdArgs);
		
		assertEquals(1.0f, filteredSNPs.getScore(/*traitIdx*/0, /*snpIdx*/0), 1e-3);
		assertEquals(1.1f, filteredSNPs.getScore(/*traitIdx*/0, /*snpIdx*/1), 1e-3);
		assertEquals(2.4f, filteredSNPs.getScore(/*traitIdx*/0, /*snpIdx*/2), 1e-3);
		assertNull(filteredSNPs.getScore(/*traitIdx*/0, /*snpIdx*/3));
		assertEquals(2.1f, filteredSNPs.getScore(/*traitIdx*/0, /*snpIdx*/4), 1e-3);
		assertEquals(1.5f, filteredSNPs.getScore(/*traitIdx*/0, /*snpIdx*/5), 1e-3);
		assertEquals(3.0f, filteredSNPs.getScore(/*traitIdx*/0, /*snpIdx*/6), 1e-3);
		assertEquals(4.1f, filteredSNPs.getScore(/*traitIdx*/0, /*snpIdx*/7), 1e-3);
		assertNull(filteredSNPs.getScore(/*traitIdx*/0, /*snpIdx*/8));
		assertEquals(0.7f, filteredSNPs.getScore(/*traitIdx*/0, /*snpIdx*/9), 1e-3);
		assertEquals(3.4f, filteredSNPs.getScore(/*traitIdx*/1, /*snpIdx*/0), 1e-3);
		assertEquals(5.8f, filteredSNPs.getScore(/*traitIdx*/1, /*snpIdx*/1), 1e-3);
		assertEquals(1.9f, filteredSNPs.getScore(/*traitIdx*/1, /*snpIdx*/2), 1e-3);
		assertEquals(4.3f, filteredSNPs.getScore(/*traitIdx*/1, /*snpIdx*/3), 1e-3);
		assertEquals(1.2f, filteredSNPs.getScore(/*traitIdx*/1, /*snpIdx*/4), 1e-3);
		assertEquals(7.2f, filteredSNPs.getScore(/*traitIdx*/1, /*snpIdx*/5), 1e-3);
		assertEquals(2.2f, filteredSNPs.getScore(/*traitIdx*/1, /*snpIdx*/6), 1e-3);
		assertEquals(3.3f, filteredSNPs.getScore(/*traitIdx*/1, /*snpIdx*/7), 1e-3);
		assertNull(filteredSNPs.getScore(/*traitIdx*/1, /*snpIdx*/8));
		assertEquals(6.4f, filteredSNPs.getScore(/*traitIdx*/1, /*snpIdx*/9), 1e-3);
		
		assertEquals(1, filteredSNPs.getNumLociNoQScore());  // SNP9 doesn't appear in the score file, so its q-score is not checked.
		assertEquals(2, filteredSNPs.getNumMonoLoci());  // SNP4 and SNP8
		assertEquals(2, filteredSNPs.getNumAmbiguousLoci());  // SNP2 and SNP7
		assertEquals(2, filteredSNPs.getNumLociNoScore(0));
		assertEquals(1, filteredSNPs.getNumLociNoScore(1));
		
		assertEquals(0, filteredSNPs.getMatchNum(AlleleMatchScheme.MATCH_NONE));
		assertEquals(3, filteredSNPs.getMatchNum(AlleleMatchScheme.MATCH_ALLELE1));  // SNP1, SNP2, SNP5
		assertEquals(0, filteredSNPs.getMatchNum(AlleleMatchScheme.MATCH_ALLELE1_FLIPPED));
		assertEquals(2, filteredSNPs.getMatchNum(AlleleMatchScheme.MATCH_ALLELE2));  // SNP3, SNP7
		assertEquals(2, filteredSNPs.getMatchNum(AlleleMatchScheme.MATCH_ALLELE2_FLIPPED));  // SNP6, SNP10
		
		assertEquals(AlleleMatchScheme.MATCH_ALLELE1, filteredSNPs.getAlleleMatchScheme(0));
		assertEquals(AlleleMatchScheme.MATCH_ALLELE1, filteredSNPs.getAlleleMatchScheme(1));
		assertEquals(AlleleMatchScheme.MATCH_ALLELE2, filteredSNPs.getAlleleMatchScheme(2));
		assertEquals(AlleleMatchScheme.MATCH_NONE, filteredSNPs.getAlleleMatchScheme(3));
		assertEquals(AlleleMatchScheme.MATCH_ALLELE1, filteredSNPs.getAlleleMatchScheme(4));
		assertEquals(AlleleMatchScheme.MATCH_ALLELE2_FLIPPED, filteredSNPs.getAlleleMatchScheme(5));
		assertEquals(AlleleMatchScheme.MATCH_ALLELE2, filteredSNPs.getAlleleMatchScheme(6));
		assertEquals(AlleleMatchScheme.MATCH_NONE, filteredSNPs.getAlleleMatchScheme(7));
		assertEquals(AlleleMatchScheme.MATCH_NONE, filteredSNPs.getAlleleMatchScheme(8));
		assertEquals(AlleleMatchScheme.MATCH_ALLELE2_FLIPPED, filteredSNPs.getAlleleMatchScheme(9));
		
		assertEquals(4, filteredSNPs.getNumLocusGroups());
		
		assertTrue(filteredSNPs.isInLocusGroup(0, 0));
		assertTrue(!filteredSNPs.isInLocusGroup(0, 1));
		assertTrue(!filteredSNPs.isInLocusGroup(0, 2));
		assertTrue(!filteredSNPs.isInLocusGroup(0, 3));
		assertTrue(!filteredSNPs.isInLocusGroup(1, 0));
		assertTrue(filteredSNPs.isInLocusGroup(1, 1));
		assertTrue(!filteredSNPs.isInLocusGroup(1, 2));
		assertTrue(!filteredSNPs.isInLocusGroup(1, 3));
		assertTrue(!filteredSNPs.isInLocusGroup(2, 0));
		assertTrue(!filteredSNPs.isInLocusGroup(2, 1));
		assertTrue(filteredSNPs.isInLocusGroup(2, 2));
		assertTrue(!filteredSNPs.isInLocusGroup(2, 3));
		assertTrue(!filteredSNPs.isInLocusGroup(3, 0));  // SNP4 is monolithic
		assertTrue(!filteredSNPs.isInLocusGroup(3, 1));
		assertTrue(!filteredSNPs.isInLocusGroup(3, 2));
		assertTrue(!filteredSNPs.isInLocusGroup(3, 3));
		assertTrue(!filteredSNPs.isInLocusGroup(4, 0));
		assertTrue(!filteredSNPs.isInLocusGroup(4, 1));
		assertTrue(!filteredSNPs.isInLocusGroup(4, 2));
		assertTrue(filteredSNPs.isInLocusGroup(4, 3));
		assertTrue(!filteredSNPs.isInLocusGroup(5, 0));
		assertTrue(!filteredSNPs.isInLocusGroup(5, 1));
		assertTrue(!filteredSNPs.isInLocusGroup(5, 2));
		assertTrue(filteredSNPs.isInLocusGroup(5, 3));
		assertTrue(!filteredSNPs.isInLocusGroup(6, 0));
		assertTrue(!filteredSNPs.isInLocusGroup(6, 1));
		assertTrue(!filteredSNPs.isInLocusGroup(6, 2));
		assertTrue(!filteredSNPs.isInLocusGroup(6, 3));
		assertTrue(!filteredSNPs.isInLocusGroup(7, 0));
		assertTrue(!filteredSNPs.isInLocusGroup(7, 1));
		assertTrue(!filteredSNPs.isInLocusGroup(7, 2));
		assertTrue(!filteredSNPs.isInLocusGroup(7, 3));
		assertTrue(!filteredSNPs.isInLocusGroup(8, 0));
		assertTrue(!filteredSNPs.isInLocusGroup(8, 1));
		assertTrue(!filteredSNPs.isInLocusGroup(8, 2));
		assertTrue(!filteredSNPs.isInLocusGroup(8, 3));
		assertTrue(!filteredSNPs.isInLocusGroup(9, 0));
		assertTrue(filteredSNPs.isInLocusGroup(9, 1));
		assertTrue(filteredSNPs.isInLocusGroup(9, 2));
		assertTrue(!filteredSNPs.isInLocusGroup(9, 3));
		
		assertEquals(1, filteredSNPs.getNumInLocusGroup(0));
		assertEquals(2, filteredSNPs.getNumInLocusGroup(1));
		assertEquals(2, filteredSNPs.getNumInLocusGroup(2));
		assertEquals(2, filteredSNPs.getNumInLocusGroup(3));
	}

}
