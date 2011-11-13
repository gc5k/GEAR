package family.pedigree.design.hierarchy;

import java.util.ArrayList;
import java.util.Hashtable;

import org.apache.commons.lang3.ArrayUtils;

import util.NewIt;
import util.Sample;
import family.mdr.data.PersonIndex;
import family.pedigree.file.PhenotypeFile;
import family.pedigree.file.MapFile;
import family.pedigree.file.PedigreeFile;
import family.pedigree.genotype.BFamilyStruct;
import family.pedigree.genotype.BPerson;
import family.pedigree.phenotype.FamilyUnit;
import family.pedigree.phenotype.Subject;

/**
 * 
 * @author Guo-Bo Chen, chenguobo@gmail.com
 */
public final class SII extends ChenBase {

	public SII(PedigreeFile ped, PhenotypeFile phe, MapFile map, long s, int pIdx, int[] cIdx, int m) {
		super(ped, phe, map, s, pIdx, cIdx, m);
	}

	public void RevvingUp() {

		ChenBase.LineUpGenotypePhenotype lineup = new ChenBase.LineUpGenotypePhenotype();
		Hashtable<String, BFamilyStruct> Fam = PedData.getFamilyStruct();

		PersonTable.ensureCapacity(qualified_Sib);

		if (PhenoData != null)
			CovariateTable.ensureCapacity(qualified_Sib);

		status = new double[qualified_Sib];

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
			BFamilyStruct fs = Fam.get(fi);
			FamilyUnit FamUnit = PhenoData == null ? null : PhenoData.getFamilyUnit(fi);
			String[] pi = fs.getPersonListSorted();
			int si = 0;

			for (int i = 0; i < pi.length; i++) {
				if (!lineup.filter[c][i])
					continue;
				BPerson per = fs.getPerson(pi[i]);
				Subject sub = PhenoData == null ? null : FamUnit.getSubject(pi[i]);

				if (fs.hasAncestor(per)) {
					s_P.add(new PersonIndex(fs.getFamilyStructName(), pi[i], per, false, false));
					status[s] = Double.parseDouble(per.getAffectedStatus());
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
	}

	public void getPermutedScore(boolean isNested) {

		if (isNested) {
			int c = 0;
			for (int i = 0; i < numSib.length; i++) {
				int[] si = Sample.SampleIndex(0, numSib[i] - 1, numSib[i]);
				for (int j = 0; j < si.length; j++) {
					PersonTable.get(c+j).setPermutedScore(score[c+si[j]]);
				}
				c += si.length;
			}
		} else {
			int[] idx = Sample.SampleIndex(0, PersonTable.size() - 1, PersonTable.size());
			for (int i = 0; i < idx.length; i++) {
				PersonTable.get(i).setPermutedScore(score[idx[i]]);
			}
		}
	}

}
