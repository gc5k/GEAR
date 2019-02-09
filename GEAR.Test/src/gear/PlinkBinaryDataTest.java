package gear;

import gear.data.Person;
import gear.data.Family;
import gear.data.UniqueRecordSet;
import gear.family.pedigree.file.MapFile;
import gear.family.pedigree.file.PedigreeFile;
import gear.family.pedigree.file.SNP;
import gear.family.plink.PLINKBinaryParser;

import static org.junit.Assert.*;

import org.junit.Test;

public class PlinkBinaryDataTest
{
	@Test
	public void testSnpMajor()
	{
		test("data/sim.bed", "data/sim.bim", "data/sim.fam");
	}
	
	@Test
	public void testIndividualMajor()
	{
		test("data/sim_ind_major.bed", "data/sim.bim", "data/sim.fam");
	}
	
	public void test(String bed, String bim, String fam)
	{
		PLINKBinaryParser parser = new PLINKBinaryParser(bed, bim, fam);
		parser.parse();
		
		PedigreeFile pedData = parser.getPedigreeData();
		
		assertEquals(6, pedData.getNumIndividuals());
		assertEquals(5, pedData.getNumMarkers());

		// individual info begin
		UniqueRecordSet<Family> families = pedData.getFamilies();
		assertEquals(4, families.size());
		
		Family fam0 = families.get("fam0"); 
		assertNotNull(fam0);
		assertSame(fam0, families.get(0));
		assertEquals(3, fam0.size());
		
		Person per0 = fam0.getPerson("per0");
		assertNotNull(per0);
		assertSame(per0, fam0.getPerson(0));
		assertEquals("0", per0.getDadID());
		assertEquals("0", per0.getMomID());
		assertEquals(1, per0.getGender());
		assertEquals("2", per0.getAffectedStatus());
		
		Person per1 = fam0.getPerson("per1");
		assertNotNull(per1);
		assertSame(per1, fam0.getPerson(1));
		assertEquals("0", per1.getDadID());
		assertEquals("0", per1.getMomID());
		assertEquals(2, per1.getGender());
		assertEquals("2", per1.getAffectedStatus());
		
		Person per2 = fam0.getPerson("per2");
		assertNotNull(per2);
		assertSame(per2, fam0.getPerson(2));
		assertEquals("per0", per2.getDadID());
		assertEquals("per1", per2.getMomID());
		assertEquals(2, per2.getGender());
		assertEquals("2", per2.getAffectedStatus());
		
		Family per3Fam = families.get("per3"); 
		assertNotNull(per3Fam);
		assertSame(per3Fam, families.get(1));
		assertEquals(1, per3Fam.size());
		
		Person per3 = per3Fam.getPerson("per3");
		assertNotNull(per3);
		assertSame(per3, per3Fam.getPerson(0));
		assertEquals("0", per3.getDadID());
		assertEquals("0", per3.getMomID());
		assertEquals(2, per3.getGender());
		assertEquals("2", per3.getAffectedStatus());
		
		Family per4Fam = families.get("per4"); 
		assertNotNull(per4Fam);
		assertSame(per4Fam, families.get(2));
		assertEquals(1, per4Fam.size());
		
		Person per4 = per4Fam.getPerson("per4");
		assertNotNull(per4);
		assertSame(per4, per4Fam.getPerson(0));
		assertEquals("0", per4.getDadID());
		assertEquals("0", per4.getMomID());
		assertEquals(2, per4.getGender());
		assertEquals("1", per4.getAffectedStatus());
		
		Family per5Fam = families.get("per5"); 
		assertNotNull(per5Fam);
		assertSame(per5Fam, families.get(3));
		assertEquals(1, per5Fam.size());
		
		Person per5 = per5Fam.getPerson("per5");
		assertNotNull(per5);
		assertSame(per5, per5Fam.getPerson(0));
		assertEquals("0", per5.getDadID());
		assertEquals("0", per5.getMomID());
		assertEquals(2, per5.getGender());
		assertEquals("1", per5.getAffectedStatus());
		// individual info end
		
		// locus info begin
		MapFile mapData = parser.getMapData();
		assertEquals(5, mapData.getMarkerNumber());
		
		SNP snp0 = mapData.getSNP(0);
		assertEquals("1", snp0.getChromosome());
		assertEquals("null_0", snp0.getName());
		assertEquals(1, snp0.getPosition());
		assertEquals('D', snp0.getSNP()[0]);
		assertEquals('d', snp0.getSNP()[1]);
		
		SNP snp1 = mapData.getSNP(1);
		assertEquals("1", snp1.getChromosome());
		assertEquals("null_1", snp1.getName());
		assertEquals(2, snp1.getPosition());
		assertEquals('d', snp1.getSNP()[0]);
		assertEquals('D', snp1.getSNP()[1]);
		
		SNP snp2 = mapData.getSNP(2);
		assertEquals("1", snp2.getChromosome());
		assertEquals("null_2", snp2.getName());
		assertEquals(3, snp2.getPosition());
		assertEquals('d', snp2.getSNP()[0]);
		assertEquals('D', snp2.getSNP()[1]);
		
		SNP snp3 = mapData.getSNP(3);
		assertEquals("1", snp3.getChromosome());
		assertEquals("null_3", snp3.getName());
		assertEquals(4, snp3.getPosition());
		assertEquals('D', snp3.getSNP()[0]);
		assertEquals('d', snp3.getSNP()[1]);
		
		SNP snp4 = mapData.getSNP(4);
		assertEquals("1", snp4.getChromosome());
		assertEquals("null_4", snp4.getName());
		assertEquals(5, snp4.getPosition());
		assertEquals('d', snp4.getSNP()[0]);
		assertEquals('D', snp4.getSNP()[1]);
		// locus info end
		
		// genotypes begin
		assertEquals("22", per0.getBiAlleleGenotypeString(0));
		assertEquals("01", per0.getBiAlleleGenotypeString(1));
		assertEquals("01", per0.getBiAlleleGenotypeString(2));
		assertEquals("01", per0.getBiAlleleGenotypeString(3));
		assertEquals("10", per0.getBiAlleleGenotypeString(4));
		
		assertEquals("10", per1.getBiAlleleGenotypeString(0));
		assertEquals("10", per1.getBiAlleleGenotypeString(1));
		assertEquals("10", per1.getBiAlleleGenotypeString(2));
		assertEquals("01", per1.getBiAlleleGenotypeString(3));
		assertEquals("10", per1.getBiAlleleGenotypeString(4));
		
		assertEquals("10", per2.getBiAlleleGenotypeString(0));
		assertEquals("10", per2.getBiAlleleGenotypeString(1));
		assertEquals("01", per2.getBiAlleleGenotypeString(2));
		assertEquals("01", per2.getBiAlleleGenotypeString(3));
		assertEquals("01", per2.getBiAlleleGenotypeString(4));
		
		assertEquals("01", per3.getBiAlleleGenotypeString(0));
		assertEquals("10", per3.getBiAlleleGenotypeString(1));
		assertEquals("01", per3.getBiAlleleGenotypeString(2));
		assertEquals("10", per3.getBiAlleleGenotypeString(3));
		assertEquals("01", per3.getBiAlleleGenotypeString(4));
		
		assertEquals("00", per4.getBiAlleleGenotypeString(0));
		assertEquals("10", per4.getBiAlleleGenotypeString(1));
		assertEquals("01", per4.getBiAlleleGenotypeString(2));
		assertEquals("01", per4.getBiAlleleGenotypeString(3));
		assertEquals("10", per4.getBiAlleleGenotypeString(4));
		
		assertEquals("10", per5.getBiAlleleGenotypeString(0));
		assertEquals("01", per5.getBiAlleleGenotypeString(1));
		assertEquals("01", per5.getBiAlleleGenotypeString(2));
		assertEquals("01", per5.getBiAlleleGenotypeString(3));
		assertEquals("10", per5.getBiAlleleGenotypeString(4));
		// genotypes end
	}
}
