package gear.family.pedigree.file;

import static org.junit.Assert.*;

import org.junit.Before;
import org.junit.Test;

public class BEDReaderTest {
	
	private BEDReader bedReader;

	@Before
	public void setUp() throws Exception {
		BIMReader bim = new BIMReader("data/sim.bim");
		bedReader = new BEDReader("data/sim.bed", "data/sim.fam", bim);
		bedReader.prepareToParseGenotypes();
	}

	@Test
	public void test() {
		assertEquals(0b10111101, bedReader.readNextByte());
		bedReader.reset();
		assertEquals(0b10111101, bedReader.readNextByte());
		assertEquals(0b00001100, bedReader.readNextByte());
		bedReader.skipOneRow();
		assertEquals(0b10101110, bedReader.readNextByte());
		assertEquals(0b00001010, bedReader.readNextByte());
	}

}
