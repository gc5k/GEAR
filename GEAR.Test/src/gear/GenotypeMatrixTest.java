package gear;

import static org.junit.Assert.*;

import java.util.ArrayList;

import org.junit.Test;

import gear.data.Person;
import gear.family.pedigree.PersonIndex;
import gear.family.popstat.GenotypeMatrix;

public class GenotypeMatrixTest
{

	@Test
	public void testAdditiveScoreGetterSetter()
	{
		Person person = new Person(31 /*markers*/);
		PersonIndex personIdx = new PersonIndex("FAM1", "IND1", person, false /*pseudo*/, false /*founder*/);
		ArrayList<PersonIndex> personIndexes = new ArrayList<PersonIndex>();
		personIndexes.add(personIdx);
		GenotypeMatrix genoMat = new GenotypeMatrix(personIndexes);
		
		genoMat.setAdditiveScore(0 /*person*/, 0 /*marker*/, 0 /*genotype*/);
		genoMat.setAdditiveScore(0 /*person*/, 1 /*marker*/, 1 /*genotype*/);
		genoMat.setAdditiveScore(0 /*person*/, 2 /*marker*/, 2 /*genotype*/);
		genoMat.setAdditiveScore(0 /*person*/, 3 /*marker*/, 3 /*genotype*/);
		genoMat.setAdditiveScore(0 /*person*/, 4 /*marker*/, 2 /*genotype*/);
		genoMat.setAdditiveScore(0 /*person*/, 30 /*marker*/, 1 /*genotype*/);
		genoMat.setAdditiveScore(0 /*person*/, 31 /*marker*/, 2 /*genotype*/);
		
		assertEquals(0, genoMat.getAdditiveScore(0 /*person*/, 0 /*marker*/));
		assertEquals(1, genoMat.getAdditiveScore(0 /*person*/, 1 /*marker*/));
		assertEquals(2, genoMat.getAdditiveScore(0 /*person*/, 2 /*marker*/));
		assertEquals(3, genoMat.getAdditiveScore(0 /*person*/, 3 /*marker*/));
		assertEquals(2, genoMat.getAdditiveScore(0 /*person*/, 4 /*marker*/));
		assertEquals(1, genoMat.getAdditiveScore(0 /*person*/, 30 /*marker*/));
		assertEquals(2, genoMat.getAdditiveScore(0 /*person*/, 31 /*marker*/));
	}

}
