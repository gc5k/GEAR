package gear.subcommands.profile;

import static org.junit.Assert.*;
import gear.subcommands.profile.Score;
import gear.subcommands.profile.ScoreFile;

import org.junit.Test;

public class ScoreFileTest
{
	@Test
	public void testHasHeader()
	{
		ScoreFile scoreFile = new ScoreFile("data/ScoresWithHeaders.txt", /*hasHeaders*/true);
		
		assertEquals(2, scoreFile.getNumberOfTraits());
		assertEquals("Trait1", scoreFile.getTrait(0));
		assertEquals("Trait2", scoreFile.getTrait(1));
		
		Score score;
		
		score = scoreFile.getScore("SNP1");
		assertNotNull(score);
		assertEquals('A', score.getAllele());
		assertEquals(1.0f, score.getValue(0), 1e-3);
		assertEquals(3.4f, score.getValue(1), 1e-3);
		
		score = scoreFile.getScore("SNP2");
		assertNotNull(score);
		assertEquals('A', score.getAllele());
		assertEquals(1.1f, score.getValue(0), 1e-3);
		assertEquals(5.8f, score.getValue(1), 1e-3);
		
		score = scoreFile.getScore("SNP3");
		assertNotNull(score);
		assertEquals('G', score.getAllele());
		assertEquals(2.4f, score.getValue(0), 1e-3);
		assertEquals(1.9f, score.getValue(1), 1e-3);
		
		score = scoreFile.getScore("SNP4");
		assertNotNull(score);
		assertEquals('C', score.getAllele());
		assertNull(score.getValue(0));
		assertEquals(4.3f, score.getValue(1), 1e-3);
		
		score = scoreFile.getScore("SNP5");
		assertNotNull(score);
		assertEquals('C', score.getAllele());
		assertEquals(2.1f, score.getValue(0), 1e-3);
		assertEquals(1.2f, score.getValue(1), 1e-3);
		
		score = scoreFile.getScore("SNP6");
		assertNotNull(score);
		assertEquals('A', score.getAllele());
		assertEquals(1.5f, score.getValue(0), 1e-3);
		assertEquals(7.2f, score.getValue(1), 1e-3);
		
		score = scoreFile.getScore("SNP7");
		assertNotNull(score);
		assertEquals('G', score.getAllele());
		assertEquals(3.0f, score.getValue(0), 1e-3);
		assertEquals(2.2f, score.getValue(1), 1e-3);
		
		score = scoreFile.getScore("SNP8");
		assertNotNull(score);
		assertEquals('C', score.getAllele());
		assertEquals(4.1f, score.getValue(0), 1e-3);
		assertEquals(3.3f, score.getValue(1), 1e-3);
		
		assertNull(scoreFile.getScore("SNP9"));
		
		score = scoreFile.getScore("SNP10");
		assertNotNull(score);
		assertEquals('C', score.getAllele());
		assertEquals(0.7f, score.getValue(0), 1e-3);
		assertEquals(6.4f, score.getValue(1), 1e-3);
	}
	
	@Test
	public void testNoHeader()
	{
		ScoreFile scoreFile = new ScoreFile("data/ScoresWithoutHeaders.txt", /*hasHeaders*/false);
		
		assertEquals(2, scoreFile.getNumberOfTraits());
		assertEquals("1", scoreFile.getTrait(0));
		assertEquals("2", scoreFile.getTrait(1));
		
		Score score;
		
		score = scoreFile.getScore("SNP1");
		assertNotNull(score);
		assertEquals('A', score.getAllele());
		assertEquals(1.0f, score.getValue(0), 1e-3);
		assertEquals(3.4f, score.getValue(1), 1e-3);
		
		score = scoreFile.getScore("SNP2");
		assertNotNull(score);
		assertEquals('A', score.getAllele());
		assertEquals(1.1f, score.getValue(0), 1e-3);
		assertEquals(5.8f, score.getValue(1), 1e-3);
		
		score = scoreFile.getScore("SNP3");
		assertNotNull(score);
		assertEquals('G', score.getAllele());
		assertEquals(2.4f, score.getValue(0), 1e-3);
		assertEquals(1.9f, score.getValue(1), 1e-3);
		
		score = scoreFile.getScore("SNP4");
		assertNotNull(score);
		assertEquals('C', score.getAllele());
		assertNull(score.getValue(0));
		assertEquals(4.3f, score.getValue(1), 1e-3);
		
		score = scoreFile.getScore("SNP5");
		assertNotNull(score);
		assertEquals('C', score.getAllele());
		assertEquals(2.1f, score.getValue(0), 1e-3);
		assertEquals(1.2f, score.getValue(1), 1e-3);
		
		score = scoreFile.getScore("SNP6");
		assertNotNull(score);
		assertEquals('A', score.getAllele());
		assertEquals(1.5f, score.getValue(0), 1e-3);
		assertEquals(7.2f, score.getValue(1), 1e-3);
		
		score = scoreFile.getScore("SNP7");
		assertNotNull(score);
		assertEquals('G', score.getAllele());
		assertEquals(3.0f, score.getValue(0), 1e-3);
		assertEquals(2.2f, score.getValue(1), 1e-3);
		
		score = scoreFile.getScore("SNP8");
		assertNotNull(score);
		assertEquals('C', score.getAllele());
		assertEquals(4.1f, score.getValue(0), 1e-3);
		assertEquals(3.3f, score.getValue(1), 1e-3);
		
		assertNull(scoreFile.getScore("SNP9"));
		
		score = scoreFile.getScore("SNP10");
		assertNotNull(score);
		assertEquals('C', score.getAllele());
		assertEquals(0.7f, score.getValue(0), 1e-3);
		assertEquals(6.4f, score.getValue(1), 1e-3);
	}
}
