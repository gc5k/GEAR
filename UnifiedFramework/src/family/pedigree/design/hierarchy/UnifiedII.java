package family.pedigree.design.hierarchy;

import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Random;
import java.util.Map.Entry;

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
	
	public UnifiedII(String ped, String phe) {
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
			numSibs += fs.getNumSibs();
		}

		PersonTable.ensureCapacity(numUnrelated + numSibs);
		CovariateTable.ensureCapacity(numUnrelated + numSibs);
		genotype = new byte[numUnrelated + numSibs][];
		status = new byte[numUnrelated + numSibs];

		int n_unrelated = 0;
		int n_sib = 0;

		ArrayList<PersonIndex> u_P = NewIt.newArrayList();
		ArrayList<ArrayList<String>> u_C = NewIt.newArrayList();

		ArrayList<PersonIndex> s_P = NewIt.newArrayList();
		ArrayList<ArrayList<String>> s_C = NewIt.newArrayList();

		ArrayList<Integer> SibIdx = NewIt.newArrayList();
		ArrayList<ArrayList<Integer>> FounderIdx = NewIt.newArrayList();
		ArrayList<Integer> Singleton = NewIt.newArrayList();
		
		int index = 0;
		for (String fi : PedData.getFamListSorted()) {
			FamilyStruct fs = Fam.get(fi);
			FamilyUnit FamUnit = PhenoData.getFamilyUnit(fi);
			String[] pi = fs.getPersonListSorted();
			int len = pi.length;
			int si = 0;
			ArrayList<Integer> F = NewIt.newArrayList();

			for (int i = 0; i < pi.length; i++) {
				Person per = fs.getPerson(pi[i]);
				Subject sub = FamUnit.getSubject(pi[i]);
				if (fs.hasAncestor(per) && len > 1) {
					si++;
					s_P.add(new PersonIndex(fs.getFamilyStructName(), pi[i]));
					genotype[n_sib + numUnrelated] = per.getGenotypeScore();
					status[n_sib + numUnrelated] = (byte) per.getAffectedStatus();
					s_C.add(sub.getTraits());
					n_sib++;
				} else {
					u_P.add(new PersonIndex(fs.getFamilyStructName(), pi[i]));
					genotype[n_unrelated] = per.getGenotypeScore();
					status[n_unrelated] = (byte) per.getAffectedStatus();
					u_C.add(sub.getTraits());
					n_unrelated++;
					F.add(new Integer(index++));
				}
			}
			if (si != 0)
				SibIdx.add(new Integer(si));
			if (F.size() != 0) {
				if(F.size() == 1) {
					Singleton.add(F.get(0));
				} else {
					FounderIdx.add(F);
				}
			}
		}
		PersonTable.addAll(u_P);
		PersonTable.addAll(s_P);
		CovariateTable.addAll(u_C);
		CovariateTable.addAll(s_C);

		numSib = ArrayUtils.toPrimitive(SibIdx.toArray(new Integer[0]));
		singleton = ArrayUtils.toPrimitive(Singleton.toArray(new Integer[0]));
		founder = new int[FounderIdx.size()][];
		for(int i = 0; i < founder.length; i++) {
			founder[i] = ArrayUtils.toPrimitive(FounderIdx.get(i).toArray(new Integer[0]));
		}

		AbstractGenoDistribution.rnd = new Random(seed);
		RLDriver RLD = new RLDriver();
		RLD.TDT(Fam, PedData.getMarkerInformation(), m);
	}

	public double[] getPermutedScore(boolean isNested) {
		permuted_score = new double[score.length];
		if (isNested) {
			int[] un_related = Sample.sample(singleton);//singleton
			for (int i = 0; i < un_related.length; i++) {
				permuted_score[singleton[i]] = score[un_related[i]];
			}
			
			for (int i = 0; i < founder.length; i++) {//founders
				int[] f = Sample.sample(founder[i]);
				for (int j = 0; j < f.length; j++) {
					permuted_score[founder[i][j]] = score[f[j]];
				}
			}
			
			int c = numUnrelated;
			for (int i = 0; i < numSib.length; i++) {//sibs;
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
