package gear.util;

import static org.junit.Assert.*;

import org.junit.Test;

public class BinaryInputFileTest
{

	@Test
	public void test()
	{
		BinaryInputFile file = new BinaryInputFile("data/LittleEndianTest.dat", "data");
		file.setLittleEndian(true);
		assertEquals(12, file.available());
		assertEquals(1.8f, file.readFloat(), 1e-6);
		assertEquals(8, file.available());
		assertEquals(-4.3f, file.readFloat(), 1e-6);
		assertEquals(4, file.available());
		assertEquals(2.0f, file.readFloat(), 1e-6);
		assertEquals(0, file.available());
		file.close();
	}

}
