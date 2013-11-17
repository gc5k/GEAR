package gear.data;

import static org.junit.Assert.*;

import org.junit.Test;

public class SubjectIDFileTest
{
	@Test
	public void test()
	{
		SubjectIDFile file = new SubjectIDFile("data/SubjectID.txt");
		
		assertEquals("data/SubjectID.txt", file.getFileName());
		assertEquals("file 'data/SubjectID.txt'", file.getSubjectOrderName());
		assertEquals(3, file.getNumberOfSubjects());
		
		SubjectID[] subjectIDs =
			{
				new SubjectID("FAM1", "IND1"),
				new SubjectID("FAM1", "IND2"),
				new SubjectID("FAM2", "IND1"),
				new SubjectID("FAM2", "IND2")
			};
		
		assertEquals(0, file.getSubjectIndex(subjectIDs[0]));
		assertEquals(2, file.getSubjectIndex(subjectIDs[1]));
		assertEquals(1, file.getSubjectIndex(subjectIDs[2]));
		assertEquals(-1, file.getSubjectIndex(subjectIDs[3]));
		
		assertEquals(subjectIDs[0], file.getSubjectID(0));
		assertEquals(subjectIDs[2], file.getSubjectID(1));
		assertEquals(subjectIDs[1], file.getSubjectID(2));
		
		file.swapSubjects(1, 2);
		
		assertEquals(0, file.getSubjectIndex(subjectIDs[0]));
		assertEquals(1, file.getSubjectIndex(subjectIDs[1]));
		assertEquals(2, file.getSubjectIndex(subjectIDs[2]));
		assertEquals(-1, file.getSubjectIndex(subjectIDs[3]));
		
		assertEquals(subjectIDs[0], file.getSubjectID(0));
		assertEquals(subjectIDs[1], file.getSubjectID(1));
		assertEquals(subjectIDs[2], file.getSubjectID(2));
	}
}
