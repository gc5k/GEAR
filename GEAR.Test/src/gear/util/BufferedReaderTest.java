package gear.util;

import static org.junit.Assert.*;

import org.junit.Test;

public class BufferedReaderTest
{
	@Test
	public void test()
	{
		BufferedReader reader = BufferedReader.openTextFile("data/Foo.txt", "foo");
		String[] tokens;
		tokens = reader.readTokens();
		assertNotNull(tokens);
		assertEquals(3, reader.getCurLineNum());
		assertEquals(3, tokens.length);
		assertEquals("1", tokens[0]);
		assertEquals("\\", tokens[1]);
		assertEquals(".4", tokens[2]);
		tokens = reader.readTokens();
		assertNull(tokens);
		reader.close();
	}
}
