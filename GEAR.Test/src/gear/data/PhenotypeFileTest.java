package gear.data;

import static org.junit.Assert.*;
import gear.ConstValues;

import org.junit.Test;

public class PhenotypeFileTest
{
	@Test
	public void testWithoutHeaders()
	{
		PhenotypeFile file = new PhenotypeFile("data/PhenotypeWithoutHeaders.txt", ConstValues.NO_HEADER);
		
		assertEquals("data/PhenotypeWithoutHeaders.txt", file.getFileName());
		assertEquals("file 'data/PhenotypeWithoutHeaders.txt'", file.getSubjectOrderName());
		
		testCommon(file);
	}
	
	@Test
	public void testWithHeaders()
	{
		PhenotypeFile file = new PhenotypeFile("data/PhenotypeWithHeaders.txt", /*hasHeaders*/true);
		
		assertEquals("data/PhenotypeWithHeaders.txt", file.getFileName());
		assertEquals("file 'data/PhenotypeWithHeaders.txt'", file.getSubjectOrderName());
		
		testCommon(file);
		
		assertEquals("Trait1", file.getTraitName(0));
		assertEquals("Trait2", file.getTraitName(1));
	}
	
	private void testCommon(PhenotypeFile file)
	{
		assertEquals(3, file.getNumberOfSubjects());
		assertEquals(2, file.getNumberOfTraits());
		
		SubjectID[] subjectIDs =
			{
				new SubjectID("FAM1", "IND1"),
				new SubjectID("FAM1", "IND2"),
				new SubjectID("FAM2", "IND1")
			};
		
		assertEquals(subjectIDs[0], file.getSubjectID(0));
		assertEquals(subjectIDs[1], file.getSubjectID(1));
		assertEquals(subjectIDs[2], file.getSubjectID(2));
		
		assertEquals(0, file.getSubjectIndex(subjectIDs[0]));
		assertEquals(1, file.getSubjectIndex(subjectIDs[1]));
		assertEquals(2, file.getSubjectIndex(subjectIDs[2]));
		
		assertFalse(file.isMissing(0, 0));
		assertEquals(0.4f, file.getPhenotype(0, 0), 1e-6);
		assertTrue(file.isMissing(0, 1));
		assertTrue(file.isMissing(1, 0));
		assertFalse(file.isMissing(1, 1));
		assertEquals(1.3f, file.getPhenotype(1, 1), 1e-6);
		assertFalse(file.isMissing(2, 0));
		assertEquals(-7f, file.getPhenotype(2, 0), 1e-6);
		assertFalse(file.isMissing(2, 1));
		assertEquals(10f, file.getPhenotype(2, 1), 1e-6);
		
		file.swapSubjects(0, 1);
		
		assertEquals(subjectIDs[1], file.getSubjectID(0));
		assertEquals(subjectIDs[0], file.getSubjectID(1));
		assertEquals(subjectIDs[2], file.getSubjectID(2));
		
		assertEquals(1, file.getSubjectIndex(subjectIDs[0]));
		assertEquals(0, file.getSubjectIndex(subjectIDs[1]));
		assertEquals(2, file.getSubjectIndex(subjectIDs[2]));
		
		assertTrue(file.isMissing(0, 0));
		assertFalse(file.isMissing(0, 1));
		assertEquals(1.3f, file.getPhenotype(0, 1), 1e-6);
		assertFalse(file.isMissing(1, 0));
		assertEquals(0.4f, file.getPhenotype(1, 0), 1e-6);
		assertTrue(file.isMissing(1, 1));
		assertFalse(file.isMissing(2, 0));
		assertEquals(-7f, file.getPhenotype(2, 0), 1e-6);
		assertFalse(file.isMissing(2, 1));
		assertEquals(10f, file.getPhenotype(2, 1), 1e-6);
	}
}
