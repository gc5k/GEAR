package family.pedigree.design.hierarchy;

import java.util.ArrayList;
import java.util.Hashtable;

import org.apache.commons.lang3.ArrayUtils;

import util.NewIt;
import util.Sample;
import family.mdr.data.PersonIndex;
import family.pedigree.file.GMDRPhenoFile;
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

public class UnifiedII extends ChenBase {

	private int[] singleton;
	private int[][] founder;

	public UnifiedII(PedigreeFile ped, GMDRPhenoFile phe, MapFile map, long s, int pIdx, int[] cIdx, int m) {
		super(ped, phe, map, s, pIdx, cIdx, m);
	}

	protected void RevvingUp() {

		ChenBase.LineUpGenotypePhenotype lineup = new ChenBase.LineUpGenotypePhenotype();

		Hashtable<String, BFamilyStruct> Fam = PedData.getFamilyStruct();

		PersonTable.ensureCapacity(qualified_Unrelated + qualified_Sib);
		if(PhenoData != null)
			CovariateTable.ensureCapacity(qualified_Unrelated + qualified_Sib);

//		genotype = new byte[qualified_Unrelated + qualified_Sib][];
		status = new byte[qualified_Unrelated + qualified_Sib];

		ArrayList<PersonIndex> u_P = NewIt.newArrayList();
		ArrayList<ArrayList<String>> u_C = NewIt.newArrayList();

		ArrayList<PersonIndex> s_P = NewIt.newArrayList();
		ArrayList<ArrayList<String>> s_C = NewIt.newArrayList();

		ArrayList<Integer> SibIdx = NewIt.newArrayList();
		ArrayList<ArrayList<Integer>> FounderIdx = NewIt.newArrayList();
		ArrayList<Integer> Singleton = NewIt.newArrayList();

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
			ArrayList<Integer> F = NewIt.newArrayList();

			for (int i = 0; i < pi.length; i++) {
				if (!lineup.filter[c][i])
					continue;
				BPerson per = fs.getPerson(pi[i]);
				Subject sub = PhenoData == null ? null : FamUnit.getSubject(pi[i]);

				if (fs.hasAncestor(per)) {
					s_P.add(new PersonIndex(fs.getFamilyStructName(), pi[i], per));
//					genotype[s + qualified_Unrelated] = per.getGenotypeScore();
					status[s + qualified_Unrelated] = (byte) per.getAffectedStatus();
					if(PhenoData != null)
						s_C.add(sub.getTraits());
					si++;
					s++;
				} else {
					u_P.add(new PersonIndex(fs.getFamilyStructName(), pi[i], per));
//					genotype[un] = per.getGenotypeScore();
					status[un] = (byte) per.getAffectedStatus();
					if(PhenoData != null)
						u_C.add(sub.getTraits());
					F.add(new Integer(un++));
				}
			}
			if (si != 0) {
				SibIdx.add(new Integer(si));
				if (F.size() != 0) {
					FounderIdx.add(F);
				}
			} else {
				if (F.size() == 1) {
					Singleton.add(F.get(0));
				} else {
					FounderIdx.add(F);
				}
			}
			c++;
		}
		PersonTable.addAll(u_P);
		PersonTable.addAll(s_P);
		if(PhenoData != null) {
			CovariateTable.addAll(u_C);
			CovariateTable.addAll(s_C);
		}

		numSib = ArrayUtils.toPrimitive(SibIdx.toArray(new Integer[0]));
		singleton = ArrayUtils.toPrimitive(Singleton.toArray(new Integer[0]));
		founder = new int[FounderIdx.size()][];
		for (int i = 0; i < founder.length; i++) {
			founder[i] = ArrayUtils.toPrimitive(FounderIdx.get(i).toArray(new Integer[0]));
		}
		int[] m = new int[MapData.getMarkerNumber()];
		for (int i = 0; i < m.length; i++) {
			m[i] = i;
		}
//		AbstractGenoDistribution.rnd = new Random(seed);
//		RLDriver RLD = new RLDriver();
//		RLD.TDT(Fam, getMarkerName(), m);
	}

	public void getPermutedScore(boolean isNested) {
//		permuted_score = new double[score.length];
		if (isNested) {
			int[] un_related = Sample.sample(singleton);// singleton
			for (int i = 0; i < un_related.length; i++) {
				PersonTable.get(i).setPermutedScore(score[un_related[i]]);
//				permuted_score[singleton[i]] = score[un_related[i]];
			}

			for (int i = 0; i < founder.length; i++) {// founders
				int[] f = Sample.sample(founder[i]);
				for (int j = 0; j < f.length; j++) {
					PersonTable.get(founder[i][j]).setPermutedScore(score[f[j]]);
//					permuted_score[founder[i][j]] = score[f[j]];
				}
			}

			int c = qualified_Unrelated;
			for (int i = 0; i < numSib.length; i++) {// sibs;
				int[] si = Sample.SampleIndex(0, numSib[i] - 1, numSib[i]);
				for (int j = 0; j < si.length; j++) {
					PersonTable.get(c+j).setPermutedScore(score[c + si[j]]);
//					permuted_score[c + j] = score[c + si[j]];
				}
				c += si.length;
			}
		} else {
			int[] idx = Sample.SampleIndex(0, PersonTable.size() - 1, PersonTable.size());
			for (int i = 0; i < idx.length; i++) {
				PersonTable.get(i).setPermutedScore(score[idx[i]]);
//				permuted_score[i] = score[idx[i]];
			}
		}
	}

//	@Override
//	public void RecoverScore() {
//		for (int i = 0; i < PersonTable.size(); i++) {
//			PersonTable.get(i).setPermutedScore(score[i]);
//		}
//	}
}
