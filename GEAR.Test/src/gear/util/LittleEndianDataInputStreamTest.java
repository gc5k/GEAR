package gear.util;

import static org.junit.Assert.*;

import org.junit.Test;

import java.io.*;

public class LittleEndianDataInputStreamTest
{

	@Test
	public void test()
	{
		try
		{
			LittleEndianDataInputStream strm = new LittleEndianDataInputStream(new DataInputStream(new FileInputStream("data/LittleEndianTest.dat")), Float.SIZE);
			assertTrue(strm.available() > 0);
			assertEquals(1.8f, strm.readFloat(), 1e-6);
			assertTrue(strm.available() > 0);
			assertEquals(-4.3f, strm.readFloat(), 1e-6);
			assertTrue(strm.available() > 0);
			assertEquals(2.0f, strm.readFloat(), 1e-6);
			assertEquals(strm.available(), 0);
		}
		catch (Exception e)
		{
			fail(e.getMessage());
		}
	}

}
