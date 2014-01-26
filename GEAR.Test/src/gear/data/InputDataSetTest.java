package gear.data;

import static org.junit.Assert.*;

import org.junit.Test;

public class InputDataSetTest
{
	@Test
	public void test()
	{
		InputDataSet data = new InputDataSet();
		data.readSubjectIDFile("data/SubjectID.txt");
		data.readPhenotypeFile("data/PhenotypeWithoutHeaders.txt");
		
		assertEquals(3, data.getNumberOfSubjects());
		assertEquals(2, data.getNumberOfTraits());
		
		SubjectID[] subjectIDs =
			{
				new SubjectID("FAM1", "IND1"),
				new SubjectID("FAM1", "IND2"),
				new SubjectID("FAM2", "IND1"),
				new SubjectID("FAM2", "IND2")
			};
		
		assertEquals(0, data.getSubjectIndex(subjectIDs[0]));
		assertEquals(2, data.getSubjectIndex(subjectIDs[1]));
		assertEquals(1, data.getSubjectIndex(subjectIDs[2]));
		assertEquals(-1, data.getSubjectIndex(subjectIDs[3]));
		
		assertEquals(subjectIDs[0], data.getSubjectID(0));
		assertEquals(subjectIDs[2], data.getSubjectID(1));
		assertEquals(subjectIDs[1], data.getSubjectID(2));
		
		assertFalse(data.isPhenotypeMissing(0, 0));
		assertEquals(0.4f, data.getPhenotype(0, 0), 1e-6);
		assertTrue(data.isPhenotypeMissing(0, 1));
		assertFalse(data.isPhenotypeMissing(1, 0));
		assertEquals(-7f, data.getPhenotype(1, 0), 1e-6);
		assertFalse(data.isPhenotypeMissing(1, 1));
		assertEquals(10f, data.getPhenotype(1, 1), 1e-6);
		assertTrue(data.isPhenotypeMissing(2, 0));
		assertFalse(data.isPhenotypeMissing(2, 1));
		assertEquals(1.3f, data.getPhenotype(2, 1), 1e-6);
	}
}
