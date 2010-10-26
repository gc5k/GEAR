package family.pedigree;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;

import score.CalEngine;
import score.CalEngineException;
import publicAccess.PublicData;
import edu.mit.wi.pedfile.Family;
import edu.mit.wi.pedfile.Individual;
import edu.mit.wi.pedfile.PedFileException;

public class GMDRData {
	public MDRPed PedData;
	public GMDRPhenoFile PhenoData;
	public ArrayList unMatched;
	public ArrayList Matched;

	private ArrayList RabinowitzTable; // a sampling of genotypes under null
										// distribution
	private ArrayList TransmittedTable;// Only informative Transmitted Genotypes
										// for offspring only
	private ArrayList NontransmittedTable;// corresponding Nontransmitted
											// genotypes.
	private ArrayList CovariateTable;// phenotypes of informative individuals
										// excluded parents
	private ArrayList StatusTable;// Affection status

	private ArrayList workingGenoTable;
	private ArrayList workingStatusTable;
	private ArrayList workingScoreTable;
	
	private ArrayList<Hashtable> ScoreTable;
	private ArrayList ScoreName;
	private ArrayList<PersonIndex> PersonTable;// The indexing file records the

	private Hashtable InformativePersonHash;
	private boolean IsLoadPedFile;
	private boolean IsLoadPhenoFile;
	private boolean IsGenotypeReady;

	private boolean usingFounderGenotype;
	private boolean usingChildrenGenotype;
	private boolean isLouAlgorithm;

	private ArrayList<Integer> subsetMarker;
	public static Random rnd;

	/**
	 * PersonIndex class is developed to construct a unique key for each
	 * individual.
	 * 
	 * @author Guo-Bo Chen
	 */
	public static class PersonIndex {

		String FamilyID;
		String IndividualID;
		String Key;

		/**
		 * Constructs an individual index. Key = family_id+"->"individual_in
		 * within a family.
		 * 
		 * @param fid
		 *            ID of the family where the individual comes from.
		 * @param id
		 *            ID of the individual where the individual is in a family.
		 */
		PersonIndex(String fid, String id) {
			FamilyID = new String(fid);
			IndividualID = new String(id);
			Key = new String(fid + "->" + id);
		}

		/**
		 * get the family ID of the individual
		 * 
		 * @return String
		 */
		public String getFamilyID() {
			return FamilyID;
		}

		/**
		 * get the individual ID
		 * 
		 * @return String
		 */
		public String getIndividualID() {
			return IndividualID;
		}

		/**
		 * get the key of the indivudal. Key=family_id+"->"individual_ID.
		 * 
		 * @return
		 */
		public String getKey() {
			return Key;
		}
	}

	/**
	 * Constructing the GMDR with specialized filetype
	 * 
	 * @param filetype
	 */
	public GMDRData(boolean exP, boolean usingKid, boolean isLou) {
		usingFounderGenotype = exP;
		usingChildrenGenotype = usingKid;
		isLouAlgorithm = isLou;
		PhenoData = new GMDRPhenoFile();
		PedData = new MDRPed();
		CovariateTable = new ArrayList();

		StatusTable = new ArrayList();
		ScoreTable = new ArrayList();
		ScoreName = new ArrayList();
		TransmittedTable = new ArrayList();
		PersonTable = new ArrayList();

		NontransmittedTable = new ArrayList();
		IsLoadPedFile = false;
		IsLoadPhenoFile = false;
		IsGenotypeReady = false;
		RabinowitzTable = new ArrayList();
	}

	/**
	 * Convert two alleles at a locus to a genotype value.
	 */
	public void Allele2Genotype() {
		try {
			PedData.Allele2Genotype();
			IsGenotypeReady = true;
		} catch (MDRPedFileException e) {
			e.printStackTrace(System.err);
		}
	}

	/**
	 * Impute missing genotypes
	 */
	public void GenotypeImputation() {
		try {
			PedData.GenotypeImputation();
		} catch (MDRPedFileException e) {
			e.printStackTrace(System.err);
		}
	}

	/**
	 * Build score with a selected response and predictor(s).
	 * 
	 * @param Adjust
	 *            a boolean variable. If it equals false, don't adjust the model
	 *            with predictor.
	 * @throws CalEngineException
	 */
	public void buildScore(int PIndex, int[] CIndex, boolean adjust,
			int method, boolean includeFounder, boolean includeChildren) {
		// method: 0 for regression, 1 for logistic
		Hashtable ScoreHashTable = new Hashtable();
		PersonIndex PI;
		String key;
		ArrayList<Double> T = new ArrayList();
		ArrayList<ArrayList> C = new ArrayList();
		ArrayList<PersonIndex> P = new ArrayList();

		for (int i = 0; i < PersonTable.size(); i++) {
			PI = PersonTable.get(i);
			if (!includeFounder) {
				FamilyStruct FamStr = PedData.getFamilyStruct(PI.getFamilyID());
				if (!FamStr.hasAncestor(PI.getIndividualID())) {
					continue;
				}
			}
			if (!includeChildren) {
				FamilyStruct FamStr = PedData.getFamilyStruct(PI.getFamilyID());
				if (FamStr.hasAncestor(PI.getIndividualID())) {
					continue;
				}
			}
			if (unMatched.contains(PI.getIndividualID())) {
				continue;
			}
			boolean flag = true;
			double t = 0;
			ArrayList c = new ArrayList();
			ArrayList tempc = (ArrayList) CovariateTable.get(i);
			if (PIndex == -1) {// using affecting status as phenotype
				if (((Integer) StatusTable.get(i)).intValue() == 0) {
					flag = false;
					continue;
				} else {// after converting, 0 for unaffected; 1 for affected;
						// -1 for unknown;
					t = ((Integer) StatusTable.get(i)).intValue() - 1;
				}
			} else {
				if (((String) tempc.get(PIndex))
						.compareTo(PublicData.MissingValue) == 0) {
					flag = false;
					continue;
				} else {
					t = Double.parseDouble((String) tempc.get(PIndex));
				}
			}
			if (adjust) {
				for (int j = 0; j < CIndex.length; j++) {
					if (((String) tempc.get(CIndex[j]))
							.compareTo(PublicData.MissingValue) == 0) {
						flag = false;
						break;
					}
					c.add((Double.parseDouble((String) tempc.get(CIndex[j]))));
				}
				if (!flag) {
					continue;
				}
				C.add(c);
			}
			P.add(PersonTable.get(i));
			T.add(new Double(t));
		}

		double[][] X = null;
		if (adjust) {
			X = new double[C.size()][CIndex.length];
			for (int i = 0; i < C.size(); i++) {
				ArrayList c = C.get(i);
				for (int j = 0; j < c.size(); j++) {
					X[i][j] = ((Double) c.get(j)).doubleValue();
				}
			}
		}

		double[][] Y = new double[T.size()][1];
		for (int i = 0; i < T.size(); i++) {
			Y[i][0] = ((Double) T.get(i)).doubleValue();
		}

		CalEngine CE = null;
		try {
			CE = new CalEngine(X, Y, adjust);
		} catch (CalEngineException E) {
			E.printStackTrace(System.err);
			System.exit(0);
		}
		double[][] s = CE.GeneralScore(method);
		for (int i = 0; i < P.size(); i++) {
			PI = P.get(i);
			ScoreHashTable.put(PI.getKey(), Double.toString(s[i][0]));
		}
		ScoreTable.add(ScoreHashTable);
		nameScore(PIndex, CIndex, adjust, method, includeFounder);
	}

	private void nameScore(int PIndex, int[] CIndex, boolean adjust,
			int method, boolean includeFounder) {
		String ln = new String();
		ln += method == 1 ? "linear(" : "logistic(";
		ln += PIndex == -1 ? "status->" : PhenoData.getTraitAtI(PIndex) + "->";
		for (int i = 0; i < CIndex.length; i++) {
			ln += PhenoData.getTraitAtI(CIndex[i]) + ",";
		}
		ln += includeFounder ? "included founders)" : "excluded founders)";
		ScoreName.add(ln);
	}

	/**
	 * When score has alread been built outside PedGMDR, it can be imported
	 * directly
	 * 
	 * @param score_index
	 *            Index of the score in the phenotype file
	 * @return
	 */
	public boolean fetchScore(int score_index, boolean includeFounder) {
		Hashtable ScoreHashTable = new Hashtable();
		for (int i = 0; i < PersonTable.size(); i++) {
			PersonIndex PI = PersonTable.get(i);
			FamilyStruct FamStr = PedData.getFamilyStruct(PI.getFamilyID());
			if (unMatched.contains(PI.getIndividualID())
					&& FamStr.containsPerson(PI.getIndividualID())) {
				continue;
			}
			if (score_index == -1) {
				Integer as = (Integer) StatusTable.get(i);
				if (as.intValue() == PublicData.MissingAffection) {
					continue;
				}
				ScoreHashTable.put(new String(PI.getKey()), new Double(as
						.intValue()));
			} else {
				ArrayList v = (ArrayList) CovariateTable.get(i);
				String s = (String) v.get(score_index);
				if (s.compareTo(PublicData.MissingValue) == 0) {
					continue;
				}
				ScoreHashTable.put(new String(PI.getKey()), new Double(Double
						.parseDouble(s)));
			}
		}
		if (score_index == -1) {
			ScoreName.add("status");
		} else {
			ScoreName.add(PhenoData.getTraitAtI(score_index));
		}
		return true;
	}

	private void clearTables() {
		CovariateTable.clear();
		StatusTable.clear();
		PersonTable.clear();
		TransmittedTable.clear();
		NontransmittedTable.clear();
	}

	public ArrayList getUnMatched() {
		return unMatched;
	}

	public ArrayList getMatched() {
		return Matched;
	}

	public ArrayList getRabinowitzTable() {
		return RabinowitzTable;
	}

	public ArrayList getCovariateTable() {
		return CovariateTable;
	}

	public ArrayList getStatusTable() {
		return StatusTable;
	}

	public ArrayList getNontransmittedTable() {
		return NontransmittedTable;
	}

	public ArrayList getTransmittedTable() {
		return TransmittedTable;
	}

	public ArrayList getUsedTransmittedTable(boolean isPermutation) {
		ArrayList UsedTransmittedTable = new ArrayList();
		return TransmittedTable;
	}

	public ArrayList getPopulationStatistics() {
		return PedData.getTableData();
	}

	public Hashtable getFamilyInformative() {
		return PedData.getFamInformative();
	}

	public ArrayList getMendError() {
		return PedData.getMendErrorTrace();
	}

	public ArrayList getPersonIndexTable() {
		return PersonTable;
	}

	/**
	 * Initialize basic implementation of the genotype file.
	 * 
	 * @param Ped
	 *            the name of the pedigree file
	 * @param Pheno
	 *            the name of the phenotype file
	 * @throws IOException
	 */
	public void Initial(File Ped, File Pheno) throws IOException {
		try {
			PedData.Initial(Ped);
			PedData.parseLinkage();
			PedData.check();
			PhenoData.Initial(Pheno);
			PhenoData.parsePhenotype();
		} catch (MDRPedFileException e) {
			System.err.println("PedFile initialization exception.");
			e.printStackTrace(System.err);
		} catch (GMDRPhenoFileException e) {
			System.err.println("PhenoFile initialization exception.");
			e.printStackTrace(System.err);
		} catch (PedFileException e) {
			System.err.println("Check error.");
			e.printStackTrace(System.err);
		} catch (IOException e) {
			System.err.println("GMDRData initial exception.");
			e.printStackTrace(System.err);
		}
		IsLoadPedFile = true;
		IsLoadPhenoFile = true;
	}

	public void RevvingUp(String ped, String phe) {
		ParsePedFile(ped);
		ParsePhenoFile(phe);
		Match();
	}

	/**
	 * Initialize basic implementation of the genotype file.
	 * 
	 * @param Ped
	 *            the name of the pedigree file
	 * @throws IOException
	 */
	public void ParsePedFile(String ped) {
		File PedFile = new File(ped);
		try {
			PedData.Initial(PedFile);
			PedData.parseLinkage();
		} catch (MDRPedFileException e) {
			System.err.println("PedFile initial exception.");
			e.printStackTrace(System.err);
		} catch (IOException e) {
			System.err.println("Pedgree file initialization exception.");
			e.printStackTrace(System.err);
		}
		try {
			PedData.check();
		} catch (PedFileException E) {
			E.printStackTrace(System.err);
		} catch (MDRPedFileException E) {
			E.printStackTrace(System.err);
		}
		IsLoadPedFile = true;
		Allele2Genotype();
		GenotypeImputation();
	}

	/**
	 * Initialize basic implementation of the phenotype file.
	 * 
	 * @param Pheno
	 *            the name of the pedigree file
	 * @throws IOException
	 */
	public void ParsePhenoFile(String pheno) {
		File PheFile = new File(pheno);
		try {
			PhenoData.Initial(PheFile);
			PhenoData.parsePhenotype();
		} catch (GMDRPhenoFileException e) {
			System.err.println("Pheno file initialization exception.");
			e.printStackTrace(System.err);
		} catch (IOException e) {
			System.err.println("Pheno file initialization exception.");
			e.printStackTrace(System.err);
		}
		IsLoadPhenoFile = true;
	}

	private boolean IsUnMatchedIndividual(String fid, String pid) {
		for (int i = 0; i < unMatched.size(); i++) {
			ArrayList tmp = (ArrayList) unMatched.get(i);
			if (((String) tmp.get(0)).equals(fid)
					&& ((String) tmp.get(1)).equals(pid)) {
				return true; // individual #iid in family #fid does not have
				// phenotype

			}
		}
		return false;
	}

	/**
	 * Check whether the individuals in the pedigree file match the one in the
	 * phenotype file.
	 * 
	 * @return
	 */
	public boolean Match() {
		if (IsLoadPhenoFile == false || IsLoadPedFile == false) {
			return false;
		}
		Matched = new ArrayList();
		unMatched = new ArrayList();
		int i = 0;
		Enumeration famList = PedData.getFamList();
		while (famList.hasMoreElements()) {
			i++;
			String fid = (String) famList.nextElement();
			Family Fam = PedData.getFamily(fid);
			Enumeration sibList = Fam.getMemberList();
			FamilyUnit FamUnit = PhenoData.getFamilyUnit(fid);

			while (sibList.hasMoreElements()) {
				String iid = (String) sibList.nextElement();
				try {
					Individual ind = Fam.getMember(iid);// marker

					boolean flag = false;
					Enumeration subList = FamUnit.getSubjectsList();
					while (subList.hasMoreElements()) {
						String sid = (String) subList.nextElement();
						Subject sub = FamUnit.getSubject(sid);
						if (((String) sub.getSubjectID()).equals((String) ind
								.getIndividualID())) {
							flag = true;
						}
					}
					ArrayList temp = new ArrayList();
					temp.add(fid);
					temp.add(iid);

					if (!flag) {
						System.err.println("no such individual with ID " + iid
								+ " in family " + fid
								+ " in the phenotype file.");
						unMatched.add(temp);
					} else {
						Matched.add(temp);
					}
				} catch (PedFileException e) {
					System.err.println("no such individual with ID " + iid
							+ " in family " + fid + "in the pedigree file.");
				}
			}
		}
		return true;
	}

	public void makeupWorkingTable() {
		Hashtable faminformative = PedData.getFamInformative();
		PersonIndex PI;
		int countmiss = 0;
		InformativePersonHash = new Hashtable();

		for (int i = 0; i < PersonTable.size(); i++) {
			PI = PersonTable.get(i);
			if (!((Boolean) faminformative.get(PI.getFamilyID()))
					.booleanValue()) {
				// not informative or no score
				continue;
			}
			InformativePersonHash.put(PI.getKey(), PI.IndividualID);
		}
	}

	public void RabinowitzApproach() {
		PedData.RabinowitzApproach(isLouAlgorithm);
	}

	/**
	 * It is used for real data analysis.
	 */
	public void CreateTable(boolean isRabinowitzProc) {
		if (isRabinowitzProc) {
			RabinowitzTable.clear();
		} else {
			clearTables();
		}
		if (Matched == null || unMatched == null) {
			return;
		}
		String[] FID = PedData.getFamListSorted();

		for (int i = 0; i < FID.length; i++) {
			String fid = (String) FID[i];
			FamilyStruct FamStr = PedData.getFamilyStruct(fid);
			String[] memberList = FamStr.getPersonListSorted();
			FamilyUnit FamUnit = PhenoData.getFamilyUnit(fid);
			String pid;
			Person per;
			PseudoPerson pseudoper;
			Subject sub;
			String[] PID = FamStr.getPersonListSorted();

			ArrayList<String> FounderID = new ArrayList();
			for(int j = 0; j < memberList.length; j++) {
				pid = (String) memberList[j];
				if (IsUnMatchedIndividual(fid, pid)) {
					continue;
				}
				if (!FamStr.hasAncestor(pid)) {
					FounderID.add(pid);
				}
			}
			if (isRabinowitzProc) {
				Collections.shuffle(FounderID);
			}
			int c = 0;
			for (int j = 0; j < PID.length; j++) {
				pid = PID[j];
				per = FamStr.hasAncestor(pid)?FamStr.getPerson(pid):FamStr.getPerson(FounderID.get(c++));
				pseudoper = FamStr.getPseudoPerson(pid);
				sub = FamUnit.getSubject(pid);
				if (!isRabinowitzProc) {
					PersonTable.add(new PersonIndex(fid, pid));
					CovariateTable.add((ArrayList) sub.getTraits().clone());
					StatusTable.add(new Integer(per.getAffectedStatus()));
					TransmittedTable.add(per.getGenotype().clone());
				} else {
					if (FamStr.hasAncestor(pid)) {
						RabinowitzTable.add(pseudoper.getPseudoGenotype()
								.clone());
					} else {
						RabinowitzTable.add(per.getGenotype().clone());
					}
				}
				if (isLouAlgorithm) {
					if (FamStr.hasAncestor(per.getPersonID())) {
						NontransmittedTable.add(pseudoper.getPseudoGenotype()
								.clone());
					} else {
						NontransmittedTable.add(per.getGenotype().clone());
					}
				}
			}
		}
		if (!isRabinowitzProc) {
			makeupWorkingTable();
		}
	}

	public void PrintGMDR(String TDTped, String TDTphe,
			boolean isPermutation) throws IOException {
		PrintWriter TDTpedout = new PrintWriter(TDTped);
		PrintWriter TDTpheout = isPermutation ? null : new PrintWriter(TDTphe);

//		createWorkingTable(isPermutation);
		printMissingInformation();

		ArrayList genotypeCensus = getWorkingGenoTable();
		ArrayList statusCensus = getWorkingStatusTable();
		ArrayList marker;

		ArrayList markerinfor = PedData.getMarkerInformation();
		for (int i = 0; i < markerinfor.size(); i++) {
			TDTpedout.print(markerinfor.get(i) + "\t");
		}
		TDTpedout.println("class");

		for (int i = 0; i < genotypeCensus.size(); i++) {
			marker = (ArrayList) genotypeCensus.get(i);
			for (int j = 0; j < marker.size(); j++) {
				TDTpedout.print(marker.get(j) + "\t");
			}
			TDTpedout.println(statusCensus.get(i));
		}
		TDTpedout.close();

		if (!isPermutation) {
			for (int i = 0; i < ScoreName.size(); i++) {
				if (i != ScoreName.size() - 1) {
					TDTpheout.print(ScoreName.get(i));
				} else {
					TDTpheout.print(ScoreName.get(i) + "\t");
				}
			}
			TDTpheout.println();
			ArrayList scoreConsus = getWorkingScoreTable();

			for (int j = 0; j < scoreConsus.size(); j++) {
				ArrayList sht = (ArrayList) scoreConsus.get(j);
				for (int k = 0; k < sht.size(); k++) {
					if (k != sht.size() - 1) {
						TDTpheout.print(sht.get(k) + "\t");
					} else {
						TDTpheout.print(sht.get(k));
					}
				}
				TDTpheout.println();
			}
			TDTpheout.close();
		}
	}

	public void printMissingInformation() {
		PersonIndex PI;
		int countmiss = 0;
		ArrayList Discard = new ArrayList();
		for (int i = 0; i < PersonTable.size(); i++) {
			PI = PersonTable.get(i);
			if (!InformativePersonHash.containsKey(PI.getKey())) {
				countmiss++;
				Discard.add(PI.getKey());
				continue;
			}
		}

		if (countmiss > 0) {
			System.out.println("=================");
			System.out.println("Discarded " + countmiss + " individuals.");
			for (int i = 0; i < Discard.size(); i++) {
				System.out.println((String) Discard.get(i));
			}
			System.out.println("=================");
		}
	}

	public void  CreateWorkingTable(boolean isRabinowitzProc) {
		ArrayList currentTranGeno = (isRabinowitzProc && !isLouAlgorithm) ? RabinowitzTable:TransmittedTable;
		workingGenoTable = new ArrayList();
		workingStatusTable = new ArrayList();
		workingScoreTable = new ArrayList();
		PersonIndex PI;
		String fid1 = "";
		float r = 1;
		for (int i = 0; i < PersonTable.size(); i++) {
			PI = PersonTable.get(i);
			String fid2 = PI.getFamilyID();
			FamilyStruct fs = PedData.getFamilyStruct(PI.getFamilyID());
			if (!InformativePersonHash.containsKey(PI.getKey())) {
				continue;
			}
			if (!usingChildrenGenotype && fs.hasAncestor(PI.getIndividualID())) {
				continue;
			} 
			if (!usingFounderGenotype && !fs.hasAncestor(PI.getIndividualID())) {
				continue;
			}
			if (isLouAlgorithm && fid1.compareTo(fid2) !=0) {
				r = rnd.nextFloat();
				fid1 = fid2;
			}
			int s1 = ((Integer) StatusTable.get(i)).intValue();
			String s2 = s1 == 0 ? PublicData.MissingValue:Integer.toString(s1-1);
			String s3 = s1 == 0 ? PublicData.MissingValue:Integer.toString(1-Integer.parseInt(s2));

			ArrayList tempS = new ArrayList();
			ArrayList tempSR = new ArrayList();
			for (int j =0; j < ScoreTable.size(); j++) {
				Hashtable score = (Hashtable) ScoreTable.get(j);
				if (score.containsKey(PI.getKey())) {
						tempS.add(score.get(PI.getKey()));
						double sr = -1 * Double.parseDouble((String)score.get(PI.getKey()));
						tempSR.add(Double.toString(sr));
				} else {
					tempS.add(PublicData.MissingValue);
					tempSR.add(PublicData.MissingValue);
				}
			}
			if ( r < 0.5) {
				if (isLouAlgorithm) {
					workingGenoTable.add(NontransmittedTable.get(i));
				}
				workingGenoTable.add(currentTranGeno.get(i));
				workingStatusTable.add(s2);
				workingScoreTable.add(tempS);
				if (isLouAlgorithm) {
					workingStatusTable.add(s3);
					workingScoreTable.add(tempSR);
				}
			} else {
				workingGenoTable.add(currentTranGeno.get(i));
				workingStatusTable.add(s2);
				workingScoreTable.add(tempS);
				if (isLouAlgorithm) {
					workingGenoTable.add(NontransmittedTable.get(i));
					workingStatusTable.add(s3);
					workingScoreTable.add(tempSR);
				}
			}
		}
	}

	public ArrayList getWorkingGenoTable() {
		return workingGenoTable;
	}
	
	public ArrayList getWorkingStatusTable() {
		return workingStatusTable;
	}
	
	public ArrayList<ArrayList> getWorkingScoreTable() {
		return workingScoreTable;
	}

	public ArrayList<String> getMarkerName() {
		return PedData.getMarkerInformation();
	}

	public ArrayList<String> getTraitName() {
		return PhenoData.getTraitName();
	}
	
	public void SetChosenMarker(int[] mi) {
		subsetMarker = new ArrayList();
		for(int i = 0; i < mi.length; i++) {
			subsetMarker.add(new Integer(mi[i]));
		}
	}
}
