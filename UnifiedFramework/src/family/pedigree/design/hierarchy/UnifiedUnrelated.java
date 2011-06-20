package family.pedigree.design.hierarchy;

import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Random;
import java.util.Map.Entry;

import util.NewIt;
import util.Sample;
import family.RabinowitzLairdAlgorithm.AbstractGenoDistribution;
import family.pedigree.design.RLDriver;
import family.pedigree.genotype.FamilyStruct;
import family.pedigree.genotype.Person;
import family.pedigree.phenotype.FamilyUnit;
import family.pedigree.phenotype.Subject;

public final class UnifiedUnrelated extends ChenBase {

	public UnifiedUnrelated(String ped, String phe) {
		super(ped, phe);
	}

	protected void RevvingUp(String ped, String phe) {
		ParsePedFile(ped);
		ParsePhenoFile(phe);
		int[] m = new int[PedData.getNumMarkers()];
		for (int i = 0; i < m.length; i++) {
			m[i] = i;
		}
		SetChosenMarker(m);

		Hashtable<String, FamilyStruct> Fam = PedData.getFamilyStruct();

		for (Entry<String, FamilyStruct> entry : Fam.entrySet()) {
			FamilyStruct fs = entry.getValue();
			numUnrelated += fs.getNumFounders();
		}

		PersonTable.ensureCapacity(numUnrelated);
		CovariateTable.ensureCapacity(numUnrelated);
		genotype = new byte[numUnrelated][];
		status = new byte[numUnrelated];

		int n_unrelated = 0;

		ArrayList<PersonIndex> u_P = NewIt.newArrayList();
		ArrayList<ArrayList<String>> u_C = NewIt.newArrayList();

		ArrayList<Integer> SibIdx = NewIt.newArrayList();

		for (String fi : PedData.getFamListSorted()) {
			FamilyStruct fs = Fam.get(fi);
			FamilyUnit FamUnit = PhenoData.getFamilyUnit(fi);
			String[] pi = fs.getPersonListSorted();
			int si = 0;
			for (int i = 0; i < pi.length; i++) {
				Person per = fs.getPerson(pi[i]);
				Subject sub = FamUnit.getSubject(pi[i]);
				if (!fs.hasAncestor(per)) {
					u_P.add(new PersonIndex(fs.getFamilyStructName(), pi[i]));
					genotype[n_unrelated] = per.getGenotypeScore();
					status[n_unrelated] = (byte) per.getAffectedStatus();
					u_C.add(sub.getTraits());
					n_unrelated++;
				}
			}
			if (si != 0)
				SibIdx.add(new Integer(si));
		}
		PersonTable.addAll(u_P);
		CovariateTable.addAll(u_C);

		AbstractGenoDistribution.rnd = new Random(seed);
		RLDriver RLD = new RLDriver();
		RLD.TDT(Fam, PedData.getMarkerInformation(), m);
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
