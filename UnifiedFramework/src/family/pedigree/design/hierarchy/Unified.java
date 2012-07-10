package family.pedigree.design.hierarchy;

import java.util.ArrayList;
import java.util.Hashtable;

import org.apache.commons.lang3.ArrayUtils;

import test.Test;
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

public final class Unified extends ChenBase {

	public Unified(PedigreeFile ped, PhenotypeFile phe, MapFile map, long s, int pIdx, int[] cIdx, int m) {
		super(ped, phe, map, s, pIdx, cIdx, m);
	}

	protected void RevvingUp() {

		ChenBase.LineUpGenotypePhenotype lineup = new ChenBase.LineUpGenotypePhenotype();
		Hashtable<String, BFamilyStruct> Fam = PedData.getFamilyStruct();

		PersonTable.ensureCapacity(qualified_Unrelated + qualified_Sib);
		if(PhenoData != null) 
			CovariateTable.ensureCapacity(qualified_Unrelated + qualified_Sib);

		status = new double[qualified_Unrelated + qualified_Sib];

		ArrayList<PersonIndex> u_P = NewIt.newArrayList();
		ArrayList<ArrayList<String>> u_C = NewIt.newArrayList();

		ArrayList<PersonIndex> s_P = NewIt.newArrayList();
		ArrayList<ArrayList<String>> s_C = NewIt.newArrayList();

		ArrayList<Integer> SibIdx = NewIt.newArrayList();

		int c = 0;
		int s = 0;
		int un = 0;
		for (String fi : PedData.getFamListSorted()) {
			if ((lineup.num_qualified[c][0] + lineup.num_qualified[c][1]) == 0) {
				c++;
				continue;
			}
			BFamilyStruct fs = Fam.get(fi);
			FamilyUnit FamUnit = PhenoData == null ? null:PhenoData.getFamilyUnit(fi);
			String[] pi = fs.getPersonListSorted();
			int si = 0;

			for (int i = 0; i < pi.length; i++) {
				if (!lineup.filter[c][i])
					continue;
				BPerson per = fs.getPerson(pi[i]);
				Subject sub = PhenoData == null ? null : FamUnit.getSubject(pi[i]);

				if (fs.hasAncestor(per)) {
					s_P.add(new PersonIndex(fs.getFamilyStructName(), pi[i], per, false, false));
					status[s + qualified_Unrelated] = Double.parseDouble( per.getAffectedStatus() );
					if(PhenoData != null)
						s_C.add(sub.getTraits());
					si++;
					s++;
				} else {
					u_P.add(new PersonIndex(fs.getFamilyStructName(), pi[i], per, false, true));
					status[un] = Double.parseDouble( per.getAffectedStatus() );
					if(PhenoData != null)
						u_C.add(sub.getTraits());
					un++;
				}
			}
			if (si != 0)
				SibIdx.add(new Integer(si));
			c++;
		}
		PersonTable.addAll(u_P);
		PersonTable.addAll(s_P);
		if(PhenoData != null) {
			CovariateTable.addAll(u_C);
			CovariateTable.addAll(s_C);
		}

		numSib = ArrayUtils.toPrimitive(SibIdx.toArray(new Integer[0]));

		System.err.println(un + " unrelated individuals.");
		Test.LOG.append(un + " unrelated individuals." + "\n");
		System.err.println(s + " siblings from " + numSib.length + " families.\n");
		Test.LOG.append(s + " siblings from " + numSib.length + " families." + "\n");

	}

	public void getPermutedScore(boolean isNested) {

		if (isNested) {
			int[] un_related = Sample.SampleIndex(0, qualified_Unrelated - 1, qualified_Unrelated);
			for (int i = 0; i < un_related.length; i++) {
				PersonTable.get(i).setPermutedScore(score[un_related[i]]);
			}
			int c = qualified_Unrelated;
			for (int i = 0; i < numSib.length; i++) {
				int[] si = Sample.SampleIndex(0, numSib[i] - 1, numSib[i]);
				for (int j = 0; j < si.length; j++) {
					PersonTable.get(c+j).setPermutedScore(score[c + si[j]]);
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

//	@Override
//	public void RecoverScore() {
//		for (int i = 0; i < PersonTable.size(); i++) {
//			PersonTable.get(i).setPermutedScore(score[i]);
//		}
//		
//	}
}
