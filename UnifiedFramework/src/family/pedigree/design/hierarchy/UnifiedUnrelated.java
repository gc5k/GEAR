package family.pedigree.design.hierarchy;

import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Random;

import util.NewIt;
import util.Sample;
import family.RabinowitzLairdAlgorithm.AbstractGenoDistribution;
import family.pedigree.design.RLDriver;
import family.pedigree.genotype.FamilyStruct;
import family.pedigree.genotype.Person;
import family.pedigree.phenotype.FamilyUnit;
import family.pedigree.phenotype.Subject;

/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public final class UnifiedUnrelated extends ChenBase {

	public UnifiedUnrelated(String ped, String map, String phe, long s, int pIdx, int[] cIdx, int m) {
		super(ped, map, phe, s, pIdx, cIdx, m);
	}

	protected void RevvingUp() {

		int[] m = new int[MapData.getMarkerNumber()];
		for (int i = 0; i < m.length; i++) {
			m[i] = i;
		}
		SetChosenMarker(m);
		ChenBase.LineUpGenotypePhenotype lineup = new ChenBase.LineUpGenotypePhenotype();

		Hashtable<String, FamilyStruct> Fam = PedData.getFamilyStruct();

		PersonTable.ensureCapacity(qualified_Unrelated);
		if(phenotypeFile != null)
			CovariateTable.ensureCapacity(qualified_Unrelated);

		genotype = new byte[qualified_Unrelated][];
		status = new byte[qualified_Unrelated];

		ArrayList<PersonIndex> u_P = NewIt.newArrayList();
		ArrayList<ArrayList<String>> u_C = NewIt.newArrayList();

		ArrayList<Integer> SibIdx = NewIt.newArrayList();
		int c = 0;
		int un = 0;
		for (String fi : PedData.getFamListSorted()) {
			if (lineup.num_qualified[c][0] == 0) {
				c++;
				continue;
			}
			FamilyStruct fs = Fam.get(fi);
			FamilyUnit FamUnit = phenotypeFile == null ? null:PhenoData.getFamilyUnit(fi);
			String[] pi = fs.getPersonListSorted();
			int si = 0;

			for (int i = 0; i < pi.length; i++) {
				if (!lineup.filter[c][i])
					continue;
				Person per = fs.getPerson(pi[i]);
				Subject sub = phenotypeFile == null ? null : FamUnit.getSubject(pi[i]);

				if (!fs.hasAncestor(per)) {
					u_P.add(new PersonIndex(fs.getFamilyStructName(), pi[i]));
					genotype[un] = per.getGenotypeScore();
					status[un] = (byte) per.getAffectedStatus();
					if (phenotypeFile != null)
						u_C.add(sub.getTraits());
					un++;
				}
			}
			if (si != 0)
				SibIdx.add(new Integer(si));
			c++;
		}
		PersonTable.addAll(u_P);
		if (phenotypeFile != null)
			CovariateTable.addAll(u_C);

		AbstractGenoDistribution.rnd = new Random(seed);
		RLDriver RLD = new RLDriver();
		RLD.TDT(Fam, getMarkerName(), m);
	}

	public double[] getPermutedScore(boolean isNested) {
		permuted_score = new double[score.length];
		int[] idx = Sample.SampleIndex(0, score.length - 1, score.length);
		for (int i = 0; i < idx.length; i++) {
			permuted_score[i] = score[idx[i]];
		}
		return permuted_score;
	}
}
