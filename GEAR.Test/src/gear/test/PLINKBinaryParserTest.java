package gear.test;

import java.util.*;

import junit.framework.Assert;

import family.pedigree.genotype.BFamilyStruct;
import family.pedigree.file.PedigreeFile;
import family.plink.PLINKBinaryParser;

import static org.junit.Assert.*;

import org.junit.Test;

public class PLINKBinaryParserTest {

	@Test
	public void test() {
		PLINKBinaryParser parser = new PLINKBinaryParser("data/sim.bed", "data/sim.bim", "data/sim.fam");
		parser.Parse();
		PedigreeFile pedData = parser.getPedigreeData();
		Hashtable<String, BFamilyStruct> familyStructs = pedData.getFamilyStruct();
		assertEquals(4, familyStructs.size());
		assertNotNull(familyStructs.get("fam0"));
		assertNotNull(familyStructs.get("per3"));
		assertNotNull(familyStructs.get("per4"));
		assertNotNull(familyStructs.get("per5"));
	}

}
