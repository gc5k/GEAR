package gear;

import static org.junit.Assert.*;

import org.junit.Test;

import gear.util.BufferedReader;
import gear.util.structure.MAF;

public class MAFTest
{

	@Test
	public void testNext()
	{
		BufferedReader reader = BufferedReader.openTextFile("data/maf.txt", "MAF");
		MAF maf;
		
		assertNotNull(maf = MAF.next(reader));
		assertEquals("3", maf.getChr());
		assertEquals("rs23", maf.getSNP());
		assertEquals('N', maf.getAllele1());
		assertEquals('M', maf.getAllele2());
		assertEquals(0.4f, maf.getMAF(), 1e-6);
		assertEquals(1.0f, maf.getNChr(), 1e-6);
		
		assertNotNull(maf = MAF.next(reader));
		assertEquals("7", maf.getChr());
		assertEquals("rs7", maf.getSNP());
		assertEquals('a', maf.getAllele1());
		assertEquals('b', maf.getAllele2());
		assertEquals(Double.NaN, maf.getMAF(), 1e-6);
		assertEquals(1.0f, maf.getNChr(), 1e-6);
		
		assertNull(MAF.next(reader));
	}

}
