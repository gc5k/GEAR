package gear.util;

import static org.junit.Assert.*;
import gear.data.SubjectID;

import org.junit.Test;

public class IDIndexMapTest
{
	@Test
	public void test()
	{
		IDIndexMap<SubjectID> subjectMap = new IDIndexMap<SubjectID>();
		
		SubjectID subject0 = new SubjectID("FAM1", "IND1"),
				  subject1 = new SubjectID("FAM1", "IND2"),
				  subject2 = new SubjectID("FAM2", "IND1"),
				  subject3 = new SubjectID("FAM2", "IND2");
		
		assertTrue(subjectMap.add(subject0));
		assertTrue(subjectMap.add(subject1));
		assertTrue(subjectMap.add(subject2));
		assertFalse(subjectMap.add(subject0));
		
		assertEquals(3, subjectMap.getNumberOfEntries());
		
		assertEquals(0, subjectMap.getIndex(subject0));
		assertEquals(1, subjectMap.getIndex(subject1));
		assertEquals(2, subjectMap.getIndex(subject2));
		assertEquals(-1, subjectMap.getIndex(subject3));
		
		assertEquals(subject0, subjectMap.getID(0));
		assertEquals(subject1, subjectMap.getID(1));
		assertEquals(subject2, subjectMap.getID(2));
		
		subjectMap.swapEntries(0, 2);
		assertEquals(2, subjectMap.getIndex(subject0));
		assertEquals(1, subjectMap.getIndex(subject1));
		assertEquals(0, subjectMap.getIndex(subject2));
		assertEquals(-1, subjectMap.getIndex(subject3));
		
		assertEquals(subject0, subjectMap.getID(2));
		assertEquals(subject1, subjectMap.getID(1));
		assertEquals(subject2, subjectMap.getID(0));
	}
}
