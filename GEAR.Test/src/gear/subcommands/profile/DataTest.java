package gear.subcommands.profile;

import static org.junit.Assert.*;

import org.junit.Test;

import gear.family.pedigree.file.SNP;
import gear.family.plink.PLINKBinaryParser;
import gear.subcommands.profile.MachData;

public class DataTest
{
	@Test
	public void testPlinkData()
	{
		PLINKBinaryParser parser = new PLINKBinaryParser("data/sim.bed", "data/sim.bim", "data/sim.fam");
		parser.parse();
//		gear.subcommands.profile.Data data = new gear.subcommands.profile.PlinkData(parser);
//		
//		SNP[] snps = data.getSNPs();
//		assertEquals(5, snps.length);
//		assertEquals("null_0", snps[0].getName());
//		assertEquals('D', snps[0].getFirstAllele());
//		assertEquals('d', snps[0].getSecAllele());
//		assertEquals("null_1", snps[1].getName());
//		assertEquals('d', snps[1].getFirstAllele());
//		assertEquals('D', snps[1].getSecAllele());
//		assertEquals("null_2", snps[2].getName());
//		assertEquals('d', snps[2].getFirstAllele());
//		assertEquals('D', snps[2].getSecAllele());
//		assertEquals("null_3", snps[3].getName());
//		assertEquals('D', snps[3].getFirstAllele());
//		assertEquals('d', snps[3].getSecAllele());
//		assertEquals("null_4", snps[4].getName());
//		assertEquals('d', snps[4].getFirstAllele());
//		assertEquals('D', snps[4].getSecAllele());
//		
//		gear.subcommands.profile.Data.Iterator iter = data.iterator();
//		
//		assertTrue(iter.next());
//		assertEquals(0, iter.getIndividualIndex());
//		assertEquals("fam0", iter.getFamilyID());
//		assertEquals("per0", iter.getIndividualID());
//		assertEquals("2", iter.getPhenotype());
//		assertEquals(0, iter.getLocusIndex());
//		assertEquals(2.0 * (3.0 / 10.0), iter.getAllele1Fraction(), 1e-6);
//		
//		assertTrue(iter.next());
//		assertEquals(0, iter.getIndividualIndex());
//		assertEquals("fam0", iter.getFamilyID());
//		assertEquals("per0", iter.getIndividualID());
//		assertEquals("2", iter.getPhenotype());
//		assertEquals(1, iter.getLocusIndex());
//		assertEquals(1.0, iter.getAllele1Fraction(), 1e-8);
//		
//		assertTrue(iter.next());
//		assertEquals(0, iter.getIndividualIndex());
//		assertEquals("fam0", iter.getFamilyID());
//		assertEquals("per0", iter.getIndividualID());
//		assertEquals("2", iter.getPhenotype());
//		assertEquals(2, iter.getLocusIndex());
//		assertEquals(1.0, iter.getAllele1Fraction(), 1e-8);
//		
//		assertTrue(iter.next());
//		assertEquals(0, iter.getIndividualIndex());
//		assertEquals("fam0", iter.getFamilyID());
//		assertEquals("per0", iter.getIndividualID());
//		assertEquals("2", iter.getPhenotype());
//		assertEquals(3, iter.getLocusIndex());
//		assertEquals(1.0, iter.getAllele1Fraction(), 1e-8);
//		
//		assertTrue(iter.next());
//		assertEquals(0, iter.getIndividualIndex());
//		assertEquals("fam0", iter.getFamilyID());
//		assertEquals("per0", iter.getIndividualID());
//		assertEquals("2", iter.getPhenotype());
//		assertEquals(4, iter.getLocusIndex());
//		assertEquals(0.0, iter.getAllele1Fraction(), 1e-8);
//		
//		assertTrue(iter.next());
//		assertEquals(1, iter.getIndividualIndex());
//		assertEquals("fam0", iter.getFamilyID());
//		assertEquals("per1", iter.getIndividualID());
//		assertEquals("2", iter.getPhenotype());
//		assertEquals(0, iter.getLocusIndex());
//		assertEquals(0.0, iter.getAllele1Fraction(), 1e-8);
//		
//		assertTrue(iter.next());
//		assertEquals(1, iter.getIndividualIndex());
//		assertEquals("fam0", iter.getFamilyID());
//		assertEquals("per1", iter.getIndividualID());
//		assertEquals("2", iter.getPhenotype());
//		assertEquals(1, iter.getLocusIndex());
//		assertEquals(0.0, iter.getAllele1Fraction(), 1e-8);
//		
//		assertTrue(iter.next());
//		assertEquals(1, iter.getIndividualIndex());
//		assertEquals("fam0", iter.getFamilyID());
//		assertEquals("per1", iter.getIndividualID());
//		assertEquals("2", iter.getPhenotype());
//		assertEquals(2, iter.getLocusIndex());
//		assertEquals(0.0, iter.getAllele1Fraction(), 1e-8);
//		
//		assertTrue(iter.next());
//		assertEquals(1, iter.getIndividualIndex());
//		assertEquals("fam0", iter.getFamilyID());
//		assertEquals("per1", iter.getIndividualID());
//		assertEquals("2", iter.getPhenotype());
//		assertEquals(3, iter.getLocusIndex());
//		assertEquals(1.0, iter.getAllele1Fraction(), 1e-8);
//		
//		assertTrue(iter.next());
//		assertEquals(1, iter.getIndividualIndex());
//		assertEquals("fam0", iter.getFamilyID());
//		assertEquals("per1", iter.getIndividualID());
//		assertEquals("2", iter.getPhenotype());
//		assertEquals(4, iter.getLocusIndex());
//		assertEquals(0.0, iter.getAllele1Fraction(), 1e-8);
//		
//		assertTrue(iter.next());
//		assertEquals(2, iter.getIndividualIndex());
//		assertEquals("fam0", iter.getFamilyID());
//		assertEquals("per2", iter.getIndividualID());
//		assertEquals("2", iter.getPhenotype());
//		assertEquals(0, iter.getLocusIndex());
//		assertEquals(0.0, iter.getAllele1Fraction(), 1e-8);
//		
//		assertTrue(iter.next());
//		assertEquals(2, iter.getIndividualIndex());
//		assertEquals("fam0", iter.getFamilyID());
//		assertEquals("per2", iter.getIndividualID());
//		assertEquals("2", iter.getPhenotype());
//		assertEquals(1, iter.getLocusIndex());
//		assertEquals(0.0, iter.getAllele1Fraction(), 1e-8);
//		
//		assertTrue(iter.next());
//		assertEquals(2, iter.getIndividualIndex());
//		assertEquals("fam0", iter.getFamilyID());
//		assertEquals("per2", iter.getIndividualID());
//		assertEquals("2", iter.getPhenotype());
//		assertEquals(2, iter.getLocusIndex());
//		assertEquals(1.0, iter.getAllele1Fraction(), 1e-8);
//		
//		assertTrue(iter.next());
//		assertEquals(2, iter.getIndividualIndex());
//		assertEquals("fam0", iter.getFamilyID());
//		assertEquals("per2", iter.getIndividualID());
//		assertEquals("2", iter.getPhenotype());
//		assertEquals(3, iter.getLocusIndex());
//		assertEquals(1.0, iter.getAllele1Fraction(), 1e-8);
//		
//		assertTrue(iter.next());
//		assertEquals(2, iter.getIndividualIndex());
//		assertEquals("fam0", iter.getFamilyID());
//		assertEquals("per2", iter.getIndividualID());
//		assertEquals("2", iter.getPhenotype());
//		assertEquals(4, iter.getLocusIndex());
//		assertEquals(1.0, iter.getAllele1Fraction(), 1e-8);
//		
//		assertTrue(iter.next());
//		assertEquals(3, iter.getIndividualIndex());
//		assertEquals("per3", iter.getFamilyID());
//		assertEquals("per3", iter.getIndividualID());
//		assertEquals("2", iter.getPhenotype());
//		assertEquals(0, iter.getLocusIndex());
//		assertEquals(1.0, iter.getAllele1Fraction(), 1e-8);
//		
//		assertTrue(iter.next());
//		assertEquals(3, iter.getIndividualIndex());
//		assertEquals("per3", iter.getFamilyID());
//		assertEquals("per3", iter.getIndividualID());
//		assertEquals("2", iter.getPhenotype());
//		assertEquals(1, iter.getLocusIndex());
//		assertEquals(0.0, iter.getAllele1Fraction(), 1e-8);
//		
//		assertTrue(iter.next());
//		assertEquals(3, iter.getIndividualIndex());
//		assertEquals("per3", iter.getFamilyID());
//		assertEquals("per3", iter.getIndividualID());
//		assertEquals("2", iter.getPhenotype());
//		assertEquals(2, iter.getLocusIndex());
//		assertEquals(1.0, iter.getAllele1Fraction(), 1e-8);
//		
//		assertTrue(iter.next());
//		assertEquals(3, iter.getIndividualIndex());
//		assertEquals("per3", iter.getFamilyID());
//		assertEquals("per3", iter.getIndividualID());
//		assertEquals("2", iter.getPhenotype());
//		assertEquals(3, iter.getLocusIndex());
//		assertEquals(0.0, iter.getAllele1Fraction(), 1e-8);
//		
//		assertTrue(iter.next());
//		assertEquals(3, iter.getIndividualIndex());
//		assertEquals("per3", iter.getFamilyID());
//		assertEquals("per3", iter.getIndividualID());
//		assertEquals("2", iter.getPhenotype());
//		assertEquals(4, iter.getLocusIndex());
//		assertEquals(1.0, iter.getAllele1Fraction(), 1e-8);
//		
//		assertTrue(iter.next());
//		assertEquals(4, iter.getIndividualIndex());
//		assertEquals("per4", iter.getFamilyID());
//		assertEquals("per4", iter.getIndividualID());
//		assertEquals("1", iter.getPhenotype());
//		assertEquals(0, iter.getLocusIndex());
//		assertEquals(2.0, iter.getAllele1Fraction(), 1e-8);
//		
//		assertTrue(iter.next());
//		assertEquals(4, iter.getIndividualIndex());
//		assertEquals("per4", iter.getFamilyID());
//		assertEquals("per4", iter.getIndividualID());
//		assertEquals("1", iter.getPhenotype());
//		assertEquals(1, iter.getLocusIndex());
//		assertEquals(0.0, iter.getAllele1Fraction(), 1e-8);
//		
//		assertTrue(iter.next());
//		assertEquals(4, iter.getIndividualIndex());
//		assertEquals("per4", iter.getFamilyID());
//		assertEquals("per4", iter.getIndividualID());
//		assertEquals("1", iter.getPhenotype());
//		assertEquals(2, iter.getLocusIndex());
//		assertEquals(1.0, iter.getAllele1Fraction(), 1e-8);
//		
//		assertTrue(iter.next());
//		assertEquals(4, iter.getIndividualIndex());
//		assertEquals("per4", iter.getFamilyID());
//		assertEquals("per4", iter.getIndividualID());
//		assertEquals("1", iter.getPhenotype());
//		assertEquals(3, iter.getLocusIndex());
//		assertEquals(1.0, iter.getAllele1Fraction(), 1e-8);
//		
//		assertTrue(iter.next());
//		assertEquals(4, iter.getIndividualIndex());
//		assertEquals("per4", iter.getFamilyID());
//		assertEquals("per4", iter.getIndividualID());
//		assertEquals("1", iter.getPhenotype());
//		assertEquals(4, iter.getLocusIndex());
//		assertEquals(0.0, iter.getAllele1Fraction(), 1e-8);
//		
//		assertTrue(iter.next());
//		assertEquals(5, iter.getIndividualIndex());
//		assertEquals("per5", iter.getFamilyID());
//		assertEquals("per5", iter.getIndividualID());
//		assertEquals("1", iter.getPhenotype());
//		assertEquals(0, iter.getLocusIndex());
//		assertEquals(0.0, iter.getAllele1Fraction(), 1e-8);
//		
//		assertTrue(iter.next());
//		assertEquals(5, iter.getIndividualIndex());
//		assertEquals("per5", iter.getFamilyID());
//		assertEquals("per5", iter.getIndividualID());
//		assertEquals("1", iter.getPhenotype());
//		assertEquals(1, iter.getLocusIndex());
//		assertEquals(1.0, iter.getAllele1Fraction(), 1e-8);
//		
//		assertTrue(iter.next());
//		assertEquals(5, iter.getIndividualIndex());
//		assertEquals("per5", iter.getFamilyID());
//		assertEquals("per5", iter.getIndividualID());
//		assertEquals("1", iter.getPhenotype());
//		assertEquals(2, iter.getLocusIndex());
//		assertEquals(1.0, iter.getAllele1Fraction(), 1e-8);
//		
//		assertTrue(iter.next());
//		assertEquals(5, iter.getIndividualIndex());
//		assertEquals("per5", iter.getFamilyID());
//		assertEquals("per5", iter.getIndividualID());
//		assertEquals("1", iter.getPhenotype());
//		assertEquals(3, iter.getLocusIndex());
//		assertEquals(1.0, iter.getAllele1Fraction(), 1e-8);
//		
//		assertTrue(iter.next());
//		assertEquals(5, iter.getIndividualIndex());
//		assertEquals("per5", iter.getFamilyID());
//		assertEquals("per5", iter.getIndividualID());
//		assertEquals("1", iter.getPhenotype());
//		assertEquals(4, iter.getLocusIndex());
//		assertEquals(0.0, iter.getAllele1Fraction(), 1e-8);
//		
//		assertFalse(iter.next());
//	}
//	
//	@Test
//	public void testMachData()
//	{
//		String[] dosageFiles = new String[] { "data/mach1.mldose.gz", "data/mach2.mldose.gz" };
//		String[] infoFiles = new String[] { "data/mach1.mlinfo", "data/mach2.mlinfo" };
//		
//		gear.subcommands.profile.Data data = new MachData(dosageFiles, infoFiles);
//		
//		SNP[] snps = data.getSNPs();
//		assertEquals(4, snps.length);
//		assertEquals("M1", snps[0].getName());
//		assertEquals('A', snps[0].getFirstAllele());
//		assertEquals('B', snps[0].getSecAllele());
//		assertEquals("M2", snps[1].getName());
//		assertEquals('A', snps[1].getFirstAllele());
//		assertEquals('a', snps[1].getSecAllele());
//		assertEquals("M3", snps[2].getName());
//		assertEquals('B', snps[2].getFirstAllele());
//		assertEquals('b', snps[2].getSecAllele());
//		assertEquals("M4", snps[3].getName());
//		assertEquals('1', snps[3].getFirstAllele());
//		assertEquals('2', snps[3].getSecAllele());
//		
//		gear.subcommands.profile.Data.Iterator iter = data.iterator();
//		
//		assertTrue(iter.next());
//		assertEquals(0, iter.getIndividualIndex());
//		assertEquals("FAM1", iter.getFamilyID());
//		assertEquals("IND1", iter.getIndividualID());
//		assertNull(iter.getPhenotype());
//		assertEquals(0, iter.getLocusIndex());
//		assertEquals(0.2, iter.getAllele1Fraction(), 1e-6);
//		
//		assertTrue(iter.next());
//		assertEquals(0, iter.getIndividualIndex());
//		assertEquals("FAM1", iter.getFamilyID());
//		assertEquals("IND1", iter.getIndividualID());
//		assertNull(iter.getPhenotype());
//		assertEquals(1, iter.getLocusIndex());
//		assertEquals(0.0, iter.getAllele1Fraction(), 1e-6);
//		
//		assertTrue(iter.next());
//		assertEquals(0, iter.getIndividualIndex());
//		assertEquals("FAM1", iter.getFamilyID());
//		assertEquals("IND1", iter.getIndividualID());
//		assertNull(iter.getPhenotype());
//		assertEquals(2, iter.getLocusIndex());
//		assertEquals(1.7, iter.getAllele1Fraction(), 1e-6);
//		
//		assertTrue(iter.next());
//		assertEquals(1, iter.getIndividualIndex());
//		assertEquals("FAM2", iter.getFamilyID());
//		assertEquals("IND1", iter.getIndividualID());
//		assertNull(iter.getPhenotype());
//		assertEquals(0, iter.getLocusIndex());
//		assertEquals(1.0, iter.getAllele1Fraction(), 1e-6);
//		
//		assertTrue(iter.next());
//		assertEquals(1, iter.getIndividualIndex());
//		assertEquals("FAM2", iter.getFamilyID());
//		assertEquals("IND1", iter.getIndividualID());
//		assertNull(iter.getPhenotype());
//		assertEquals(1, iter.getLocusIndex());
//		assertEquals(2.0, iter.getAllele1Fraction(), 1e-6);
//		
//		assertTrue(iter.next());
//		assertEquals(1, iter.getIndividualIndex());
//		assertEquals("FAM2", iter.getFamilyID());
//		assertEquals("IND1", iter.getIndividualID());
//		assertNull(iter.getPhenotype());
//		assertEquals(2, iter.getLocusIndex());
//		assertEquals(0.8, iter.getAllele1Fraction(), 1e-6);
//		
//		assertTrue(iter.next());
//		assertEquals(0, iter.getIndividualIndex());
//		assertEquals("FAM1", iter.getFamilyID());
//		assertEquals("IND1", iter.getIndividualID());
//		assertNull(iter.getPhenotype());
//		assertEquals(3, iter.getLocusIndex());
//		assertEquals(0.68, iter.getAllele1Fraction(), 1e-6);
//		
//		assertTrue(iter.next());
//		assertEquals(1, iter.getIndividualIndex());
//		assertEquals("FAM2", iter.getFamilyID());
//		assertEquals("IND1", iter.getIndividualID());
//		assertNull(iter.getPhenotype());
//		assertEquals(3, iter.getLocusIndex());
//		assertEquals(1.98, iter.getAllele1Fraction(), 1e-6);
//		
//		assertFalse(iter.next());
	}
}
