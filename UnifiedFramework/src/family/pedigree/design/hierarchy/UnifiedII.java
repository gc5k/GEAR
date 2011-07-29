package family.pedigree.design.hierarchy;

import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Random;

import org.apache.commons.lang3.ArrayUtils;

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

public class UnifiedII extends ChenBase {

	private int[] singleton;
	private int[][] founder;

	public UnifiedII(String ped, String map, String phe, long s, int pIdx, int[] cIdx, int m) {
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

		PersonTable.ensureCapacity(qualified_Unrelated + qualified_Sib);
		if(phenotypeFile != null)
			CovariateTable.ensureCapacity(qualified_Unrelated + qualified_Sib);

		genotype = new byte[qualified_Unrelated + qualified_Sib][];
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
			FamilyStruct fs = Fam.get(fi);
			FamilyUnit FamUnit = phenotypeFile == null ? null:PhenoData.getFamilyUnit(fi);
			String[] pi = fs.getPersonListSorted();
			int si = 0;
			ArrayList<Integer> F = NewIt.newArrayList();

			for (int i = 0; i < pi.length; i++) {
				if (!lineup.filter[c][i])
					continue;
				Person per = fs.getPerson(pi[i]);
				Subject sub = phenotypeFile == null ? null : FamUnit.getSubject(pi[i]);

				if (fs.hasAncestor(per)) {
					s_P.add(new PersonIndex(fs.getFamilyStructName(), pi[i]));
					genotype[s + qualified_Unrelated] = per.getGenotypeScore();
					status[s + qualified_Unrelated] = (byte) per.getAffectedStatus();
					if(phenotypeFile != null)
						s_C.add(sub.getTraits());
					si++;
					s++;
				} else {
					u_P.add(new PersonIndex(fs.getFamilyStructName(), pi[i]));
					genotype[un] = per.getGenotypeScore();
					status[un] = (byte) per.getAffectedStatus();
					if(phenotypeFile != null)
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
		if(phenotypeFile != null) {
			CovariateTable.addAll(u_C);
			CovariateTable.addAll(s_C);
		}

		numSib = ArrayUtils.toPrimitive(SibIdx.toArray(new Integer[0]));
		singleton = ArrayUtils.toPrimitive(Singleton.toArray(new Integer[0]));
		founder = new int[FounderIdx.size()][];
		for (int i = 0; i < founder.length; i++) {
			founder[i] = ArrayUtils.toPrimitive(FounderIdx.get(i).toArray(new Integer[0]));
		}

		AbstractGenoDistribution.rnd = new Random(seed);
		RLDriver RLD = new RLDriver();
		RLD.TDT(Fam, getMarkerName(), m);
	}

	public double[] getPermutedScore(boolean isNested) {
		permuted_score = new double[score.length];
		if (isNested) {
			int[] un_related = Sample.sample(singleton);// singleton
			for (int i = 0; i < un_related.length; i++) {
				permuted_score[singleton[i]] = score[un_related[i]];
			}

			for (int i = 0; i < founder.length; i++) {// founders
				int[] f = Sample.sample(founder[i]);
				for (int j = 0; j < f.length; j++) {
					permuted_score[founder[i][j]] = score[f[j]];
				}
			}

			int c = qualified_Unrelated;
			for (int i = 0; i < numSib.length; i++) {// sibs;
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
