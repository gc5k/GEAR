package family.pedigree.design.hierarchy;

import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Random;

import org.apache.commons.lang3.ArrayUtils;

import util.NewIt;
import util.Sample;
import family.RabinowitzLairdAlgorithm.AbstractGenoDistribution;
import family.pedigree.design.RLDriver;
import family.pedigree.file.GMDRPhenoFile;
import family.pedigree.file.MapFile;
import family.pedigree.file.PedigreeFile;
import family.pedigree.genotype.FamilyStruct;
import family.pedigree.genotype.Person;
import family.pedigree.phenotype.FamilyUnit;
import family.pedigree.phenotype.Subject;

/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public final class SII extends ChenBase {

	public SII(PedigreeFile ped, GMDRPhenoFile phe, MapFile map, long s, int pIdx, int[] cIdx, int m) {
		super(ped, phe, map, s, pIdx, cIdx, m);
	}

	public void RevvingUp() {

		ChenBase.LineUpGenotypePhenotype lineup = new ChenBase.LineUpGenotypePhenotype();
		Hashtable<String, FamilyStruct> Fam = PedData.getFamilyStruct();

		PersonTable.ensureCapacity(qualified_Sib);

		if (PhenoData != null)
			CovariateTable.ensureCapacity(qualified_Sib);

		genotype = new byte[qualified_Sib][];
		status = new byte[qualified_Sib];

		ArrayList<PersonIndex> s_P = NewIt.newArrayList();
		ArrayList<ArrayList<String>> s_C = NewIt.newArrayList();

		ArrayList<Integer> SibIdx = NewIt.newArrayList();

		int c = 0;
		int s = 0;
		for (String fi : PedData.getFamListSorted()) {
			if (lineup.num_qualified[c][1] == 0) {
				c++;
				continue;
			}
			FamilyStruct fs = Fam.get(fi);
			FamilyUnit FamUnit = PhenoData == null ? null : PhenoData.getFamilyUnit(fi);
			String[] pi = fs.getPersonListSorted();
			int si = 0;

			for (int i = 0; i < pi.length; i++) {
				if (!lineup.filter[c][i])
					continue;
				Person per = fs.getPerson(pi[i]);
				Subject sub = PhenoData == null ? null : FamUnit.getSubject(pi[i]);

				if (fs.hasAncestor(per)) {
					s_P.add(new PersonIndex(fs.getFamilyStructName(), pi[i]));
					genotype[s] = per.getGenotypeScore();
					status[s] = (byte) per.getAffectedStatus();
					if (sub != null)
						s_C.add(sub.getTraits());
					si++;
					s++;
				}
			}
			if (si != 0)
				SibIdx.add(new Integer(si));
			c++;
		}
		PersonTable.addAll(s_P);
		if (PhenoData != null)
			CovariateTable.addAll(s_C);

		numSib = ArrayUtils.toPrimitive(SibIdx.toArray(new Integer[0]));

		int[] m = new int[MapData.getMarkerNumber()];
		for (int i = 0; i < MapData.getMarkerNumber(); i++) {
			m[i] = i;
		}
		AbstractGenoDistribution.rnd = new Random(seed);
		RLDriver RLD = new RLDriver();
		RLD.TDT(Fam, getMarkerName(), m);
	}

	public double[] getPermutedScore(boolean isNested) {

		permuted_score = new double[score.length];
		if (isNested) {
			int c = 0;
			for (int i = 0; i < numSib.length; i++) {
				int[] si = Sample.SampleIndex(0, numSib[i] - 1, numSib[i]);
				for (int j = 0; j < si.length; j++) {
					permuted_score[c + j] = score[c + si[j]];
				}
				c += si.length;
			}
		} else {
			int[] idx = Sample.SampleIndex(0, score.length - 1, score.length);
			for (int i = 0; i < idx.length; i++) {
				permuted_score[i] = score[idx[i]];
			}
		}

		return permuted_score;
	}
}
