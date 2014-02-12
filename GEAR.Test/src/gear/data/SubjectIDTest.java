package gear.data;

import static org.junit.Assert.*;

import org.junit.Test;

public class SubjectIDTest
{

	@Test
	public void test()
	{
		SubjectID subject0 = new SubjectID("FAM1", "IND1");
		SubjectID subject1 = new SubjectID("FAM1", "IND1");
		assertEquals(subject0, subject1);
		assertEquals("FAM1", subject0.getFamilyID());
		assertEquals("IND1", subject0.getIndividualID());
		assertNotSame(subject0, subject1);
	}

}
