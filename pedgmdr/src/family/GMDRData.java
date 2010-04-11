package family;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Random;

import score.CalEngine;
import score.CalEngineException;
import PublicAccess.PublicData;
import edu.mit.wi.pedfile.Family;
import edu.mit.wi.pedfile.Individual;
import edu.mit.wi.pedfile.PedFileException;

public class GMDRData {

	public MDRPed PedData;
	public GMDRPhenoFile PhenoData;
	public ArrayList unMatched;
	public ArrayList Matched;
	private ArrayList TDTGenoTable;// Transmitted and Nontransmitted altogether;
	// The nontransmitted genotypes followed the corresponding transmitted
	// genotype
	private ArrayList TransmittedTable;// Only informative Transmitted Genotypes
	private ArrayList NontransmittedTable;// corresponding Nontransmitted
	// genotypes.
	private ArrayList CovariateTable;// phenotypes of informative individuals
	// excluded parents
	private ArrayList FullCovariateTable;// phenotypes of informative
	// individuals
	private ArrayList StatusTable;// Affection status
	private ArrayList<Integer> FullStatusTable;
	private ArrayList TDTScoreTable;// corresponding Score
	private ArrayList<PersonIndex> PersonTable;// The indexing file records the
	// order of individuals in
	// TransmittedTable
	// excluded parents
	private ArrayList<PersonIndex> FullPersonTable;// The indexing file records
	// the order of all
	// individuals
	private ArrayList RabinowitzTable;
	// private static String Missingvalue=".";
	private ArrayList DiscardedIndividual;
	private ArrayList CardedIndividual;
	
	private double[][] Covariates;
	private double[][] FullCovariates;
	// private double[][] Traits;
	private double[][] Status;
	private double[][] FullStatus;
	private double[][] TDTScore;
	private int[] CovIndex;
	private int[] PheIndex;
	private boolean IsLoadPedFile;
	private boolean IsLoadPhenoFile;
	private boolean PedFileType; // false for case-control design; true for
	// family based design

	private int scoreProc = 1; // 2 for logistic
	// 1 for linear

	private boolean AdjustmentOfCovariate;
	private boolean IsGenotypeReady;
	private Hashtable ScoreHashTable;

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
	public void buildScore2(int[] PIndex, int[] CIndex, boolean Adjust,
			int method, boolean includeFounder) throws CalEngineException {// method
		// 0
		// for
		// regression,
		// 1
		// for
		// logistic
		ScoreHashTable.clear();
		PersonIndex PI;
		String key;
		scoreProc = method;
		ArrayList<Double> T = new ArrayList();
		ArrayList<ArrayList> C = new ArrayList();
		ArrayList<PersonIndex> P = new ArrayList();

		setScoreProc(method);
		setAdjustmentOfCovariate(Adjust);
		setSelectedPhenoIndex(PIndex);

		if (Adjust) {
			setSelectedCovariateIndex(CIndex);
		}
		for (int i = 0; i < FullPersonTable.size(); i++) {
			PI = FullPersonTable.get(i);
			if (!includeFounder) {
				FamilyStruct FamStr = PedData.getFamilyStruct(PI.getFamilyID());
				if (!FamStr.hasAncestor(PI.getIndividualID())) {
					continue;
				}
			}
			boolean flag = true;
			double t = 0;
			ArrayList c = new ArrayList();
			ArrayList tempc = (ArrayList) FullCovariateTable.get(i);
			if (PheIndex[0] == -1) {// using affecting status as phenotype
				if (((Integer) FullStatusTable.get(i)).intValue() == 0) {// ignor
					// missing
					// value
					flag = false;
					continue;
				} else {// after converting, 0 for unaffected; 1 for affected;
					// -1 for unknown;
					t = ((Integer) FullStatusTable.get(i)).intValue() - 1;
				}
			} else {
				if (((String) tempc.get(PheIndex[0]))
						.compareTo(PublicData.MissingValue) == 0) {
					flag = false;
					continue;
				} else {
					t = Double.parseDouble((String) tempc.get(PheIndex[0]));
				}
			}
			if (Adjust) {
				for (int j = 0; j < CovIndex.length; j++) {
					if (((String) tempc.get(CovIndex[j]))
							.compareTo(PublicData.MissingValue) == 0) {
						flag = false;
						break;
					}
					c
							.add((Double.parseDouble((String) tempc
									.get(CovIndex[j]))));
				}
				if (!flag) {
					continue;
				}
				C.add(c);
			}
			P.add(FullPersonTable.get(i));
			T.add(new Double(t));
		}

		double[][] X = null;
		if (Adjust) {
			X = new double[C.size()][CovIndex.length];
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

		CalEngine CE = new CalEngine(X, Y, Adjust);
		double[][] s = CE.GeneralScore(method);
		for (int i = 0; i < P.size(); i++) {
			PI = P.get(i);
			key = PI.getKey();
			ScoreHashTable.put(new String(key), new Double(s[i][0]));
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
	public void buildScore(int[] PIndex, int[] CIndex, boolean Adjust,
			int method) throws CalEngineException {// method 0 for regression, 1
		// for logistic
		ScoreHashTable.clear();
		PersonIndex PI;
		String key;
		scoreProc = method;
		ArrayList<Double> T = new ArrayList();
		ArrayList<ArrayList> C = new ArrayList();
		ArrayList<PersonIndex> P = new ArrayList();

		setScoreProc(method);
		setAdjustmentOfCovariate(Adjust);
		setSelectedPhenoIndex(PIndex);

		if (Adjust) {
			setSelectedCovariateIndex(CIndex);
		}
		for (int i = 0; i < FullPersonTable.size(); i++) {
			boolean flag = true;
			double t = 0;
			ArrayList c = new ArrayList();
			ArrayList tempc = (ArrayList) FullCovariateTable.get(i);
			if (PheIndex[0] == -1) {// using affecting status as phenotype
				if (((Integer) FullStatusTable.get(i)).intValue() == 0) {// ignor
					// missing
					// value
					flag = false;
					continue;
				} else {// after converting, 0 for unaffected; 1 for affected;
					// -1 for unknown;
					t = ((Integer) FullStatusTable.get(i)).intValue() - 1;
				}
			} else {
				if (((String) tempc.get(PheIndex[0]))
						.compareTo(PublicData.MissingValue) == 0) {
					flag = false;
					continue;
				} else {
					t = Double.parseDouble((String) tempc.get(PheIndex[0]));
				}
			}
			if (Adjust) {
				for (int j = 0; j < CovIndex.length; j++) {
					if (((String) tempc.get(CovIndex[j]))
							.compareTo(PublicData.MissingValue) == 0) {
						flag = false;
						break;
					}
					c
							.add((Double.parseDouble((String) tempc
									.get(CovIndex[j]))));
				}
				if (!flag) {
					continue;
				}
				C.add(c);
			}
			P.add(FullPersonTable.get(i));
			T.add(new Double(t));
		}
		double[][] X = null;
		if (Adjust) {
			X = new double[C.size()][CovIndex.length];
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

		CalEngine CE = new CalEngine(X, Y, Adjust);
		double[][] s = CE.GeneralScore(method);
		for (int i = 0; i < P.size(); i++) {
			PI = P.get(i);
			key = PI.getKey();
			ScoreHashTable.put(new String(key), new Double(s[i][0]));
		}
	}

	/**
	 * Calculate score with given parameters
	 * 
	 * @param CovIndex
	 *            indexes of the covariates
	 * @param PhenoIndex
	 *            index of the response
	 * @param adjustment
	 *            whether adjust the response with covariates of not
	 * @param method
	 *            linear regression or logistic regression
	 * @throws CalEngineException
	 */
	public void CalculateScore(int[] CovIndex, int PhenoIndex,
			boolean adjustment, int method) throws CalEngineException // adjust=false
	// ,no
	// cov
	// adjust=ture, with cov
	{
		double[][] status;
		double[][] score;
		double covariate[][] = new double[FullCovariateTable.size()][CovIndex.length];

		PersonIndex PI;
		String key;
		ScoreHashTable.clear();

		for (int i = 0; i < FullCovariateTable.size(); i++) {
			for (int j = 0; j < CovIndex.length; j++) {
				covariate[i][j] = FullCovariates[i][CovIndex[j]];
			}
		}

		if (PhenoIndex == -1) {
			status = FullStatus;
		} else {
			status = new double[FullCovariates.length][1];
			for (int i = 0; i < status.length; i++) {
				for (int j = 0; j < status[i].length; j++) {
					status[i][j] = FullCovariates[i][j];
				}
			}
		}

		CalEngine CE = new CalEngine(covariate, status, adjustment);
		score = CE.GeneralScore(method);
		for (int i = 0; i < FullPersonTable.size(); i++) {
			PI = FullPersonTable.get(i);
			key = PI.getKey();
			ScoreHashTable.put(key, new Double(score[i][0]));
		}
	}

	public void CalculateScore2(int[] CovIndex, int PhenoIndex, boolean adjust,
			int method) throws CalEngineException// only sibs are used in score
	// calculation
	// adjust=false ,no cov
	// adjust=ture, with cov
	{
		double[][] status;
		double[][] score;
		double covariate[][] = new double[CovariateTable.size()][CovIndex.length];

		PersonIndex PI;
		String key;
		ScoreHashTable.clear();

		for (int i = 0; i < CovariateTable.size(); i++) {
			for (int j = 0; j < CovIndex.length; j++) {
				covariate[i][j] = Covariates[i][CovIndex[j]];
			}
		}

		if (PhenoIndex == -1) {
			status = Status;
		} else {
			status = new double[Covariates.length][1];
			for (int i = 0; i < status.length; i++) {
				for (int j = 0; j < status[i].length; j++) {
					status[i][j] = Covariates[i][j];
				}
			}
		}

		CalEngine CE = new CalEngine(covariate, status, adjust);
		score = CE.GeneralScore(method);
		for (int i = 0; i < PersonTable.size(); i++) {
			PI = PersonTable.get(i);
			key = PI.getKey();
			ScoreHashTable.put(key, new Double(score[i][0]));
		}
	}

	/**
	 * Clear all tables, including TDTGenoTable, CovariateTable,
	 * FulCovariateTable, StatusTable, TDTScoreTable, TransmittedTable,
	 * PersonTable, FullPersonTable, NontransmittedTable.
	 */
	private void clearTables() {
		TDTGenoTable.clear();
		CovariateTable.clear();
		FullCovariateTable.clear();
		StatusTable.clear();
		TDTScoreTable.clear();
		TransmittedTable.clear();
		PersonTable.clear();
		FullPersonTable.clear();
		NontransmittedTable.clear();
	}

	/**
	 * Create nontranmitted genotype with Rabinowitz & Laird algorithm.
	 */
	public void RabinowitzCreateTable() {
		RabinowitzTable.clear();
		if (Matched == null || unMatched == null) {
			return;
		}
		Enumeration famstrList = PedData.getFamStrList();
		String[] FID = new String[PedData.getNumFamilies()];
		int ind = 0;
		while (famstrList.hasMoreElements()) {
			FID[ind++] = (String) famstrList.nextElement();
		}
		Arrays.sort(FID);

		for (int i = 0; i < FID.length; i++) {
			String fid = (String) FID[i];

			FamilyStruct FamStr = PedData.getFamilyStruct(fid);

			Enumeration sibList = FamStr.getPersonList();
			FamilyUnit FamUnit = PhenoData.getFamilyUnit(fid);
			String pid;
			Person per;
			PseudoPerson pseudoper;
			Subject sub;
			while (sibList.hasMoreElements()) {
				pid = (String) sibList.nextElement();
				if (IsUnMatchedIndividual(fid, pid)) {
					continue;
				}
				try {
					per = FamStr.getPerson(pid);
					pseudoper = FamStr.getPseudoPerson(pid);
					sub = FamUnit.getSubject(pid);
					if (!FamStr.hasAncestor(per.getPersonID())) {
						continue;
					}
					if (PedFileType) {
						RabinowitzTable.add(pseudoper.getPseudoGenotype()
								.clone());
					}
				} catch (MDRPedFileException e) {
					System.err.println("Can't find the individual " + pid
							+ " genotype in family " + fid);
				} catch (GMDRPhenoFileException e) {
					System.err.println("Can't find the individual " + pid
							+ " phenotype in family " + fid);
				}
			}
		}
	}

	/**
	 * Organize data imported from the files into table forms.
	 */
	public void CreateTable() {
		if (Matched == null || unMatched == null) {
			return;
		}
		Enumeration famstrList = PedData.getFamStrList();
		String[] FID = new String[PedData.getNumFamilies()];
		int ind = 0;
		while (famstrList.hasMoreElements()) {
			FID[ind++] = (String) famstrList.nextElement();
		}
		Arrays.sort(FID);

		// while( famstrList.hasMoreElements() )
		for (int i = 0; i < FID.length; i++) {
			// String fid = (String) famstrList.nextElement();
			String fid = (String) FID[i];

			FamilyStruct FamStr = PedData.getFamilyStruct(fid);
			Enumeration sibList = FamStr.getPersonList();
			FamilyUnit FamUnit = PhenoData.getFamilyUnit(fid);
			String pid;
			Person per;
			PseudoPerson pseudoper;
			Subject sub;
			while (sibList.hasMoreElements()) {
				pid = (String) sibList.nextElement();
				if (IsUnMatchedIndividual(fid, pid)) {
					continue;
				}
				try {
					per = FamStr.getPerson(pid);
					pseudoper = FamStr.getPseudoPerson(pid);
					sub = FamUnit.getSubject(pid);
					FullPersonTable.add(new PersonIndex(fid, pid));
					FullCovariateTable.add((ArrayList) sub.getTraits().clone());
					FullStatusTable.add(new Integer(per.getAffectedStatus()));
					if (!FamStr.hasAncestor(per.getPersonID())) {
						continue;
					}

					TDTGenoTable.add(per.getGenotype().clone());
					TransmittedTable.add(per.getGenotype().clone());
					if (PedFileType) {
						TDTGenoTable.add(pseudoper.getPseudoGenotype().clone());
						NontransmittedTable.add(pseudoper.getPseudoGenotype()
								.clone());
					}
					StatusTable.add(new Integer(per.getAffectedStatus()));
					CovariateTable.add(sub.getTraits().clone());
					PersonTable.add(new PersonIndex(fid, pid));
				} catch (MDRPedFileException e) {
					System.err.println("Can't find the individual " + pid
							+ " genotype in family " + fid);
				} catch (GMDRPhenoFileException e) {
					System.err.println("Can't find the individual " + pid
							+ " phenotype in family " + fid);
				}
			}
			TableConversion();
		}
	}

	/**
	 * When score has alread been built outside PedGMDR, it can be imported
	 * directly
	 * 
	 * @param score_index
	 *            Index of the score in the phenotype file
	 * @return
	 */
	public boolean fetchScore(int score_index) {
		if (score_index < -1) {
			return false;
		} else {
			ScoreHashTable.clear();
			for (int i = 0; i < FullPersonTable.size(); i++) {
				PersonIndex pi = FullPersonTable.get(i);
				if (score_index == -1) {
					Integer as = (Integer) FullStatusTable.get(i);
					if (as.intValue() == PublicData.MissingAffection) {
						continue;
					}
					ScoreHashTable.put(new String(pi.getKey()), new Double(as
							.intValue()));
				} else {
					ArrayList v = (ArrayList) FullCovariateTable.get(i);
					String s = (String) v.get(score_index);
					if (s.compareTo(PublicData.MissingValue) == 0) {
						continue;
					}
					ScoreHashTable.put(new String(pi.getKey()), new Double(
							Double.parseDouble(s)));
				}
			}
			PheIndex = new int[1];
			PheIndex[0] = score_index;
			return true;
		}
	}

	/**
	 * Constructing the GMDR with specialized filetype
	 * 
	 * @param filetype
	 */
	public GMDRData(boolean FileType) {
		PedFileType = FileType;
		PhenoData = new GMDRPhenoFile();
		PedData = new MDRPed();
		TDTGenoTable = new ArrayList();
		CovariateTable = new ArrayList();
		FullCovariateTable = new ArrayList();
		StatusTable = new ArrayList();
		FullStatusTable = new ArrayList();
		TDTScoreTable = new ArrayList();
		TransmittedTable = new ArrayList();
		PersonTable = new ArrayList();
		FullPersonTable = new ArrayList();
		NontransmittedTable = new ArrayList();
		IsLoadPedFile = false;
		IsLoadPhenoFile = false;
		IsGenotypeReady = false;
		ScoreHashTable = new Hashtable();
		RabinowitzTable = new ArrayList();
	}

	public double[][] getCovariates() {
		return Covariates;
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

	public Hashtable getScoreHashTable() {
		return ScoreHashTable;
	}

	public ArrayList getTDTGenoTable() {
		return TDTGenoTable;
	}

	public ArrayList getCovariateTable() {
		return CovariateTable;
	}

	public ArrayList getFullCovariateTable() {
		return FullCovariateTable;
	}

	public ArrayList getStatusTable() {
		return StatusTable;
	}

	public double[][] getStatus() {
		return Status;
	}

	public ArrayList getNontransmittedTable() {
		return NontransmittedTable;
	}

	public ArrayList getTransmittedTable() {
		return TransmittedTable;
	}

	public boolean getFiletype() {
		return PedFileType;
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

	public ArrayList getFullPersonIndexTable() {
		return FullPersonTable;
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

	/**
	 * Initialize basic implementation of the genotype file.
	 * 
	 * @param Ped
	 *            the name of the pedigree file
	 * @throws IOException
	 */
	public void InitialPedFile(String ped) throws MDRPedFileException,
			PedFileException {
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
		} catch (PedFileException e) {
			System.err.println("Checking error.");
			e.printStackTrace(System.err);
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
	public void InitialPhenoFile(String pheno) throws GMDRPhenoFileException {
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
						try {
							Subject sub = FamUnit.getSubject(sid);
							if (((String) sub.getSubjectID())
									.equals((String) ind.getIndividualID())) {
								flag = true;
							}
						} catch (GMDRPhenoFileException e) {
							System.err.println("no such individual with ID "
									+ iid + " in family " + fid
									+ " in the phenotype file.");
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

	public void NonTransmittedGenoType() {
		if (PedFileType && IsGenotypeReady) {
			PedData.NonTransmittedGenoType();
		}
	}

	public double[][] pickupScore(ArrayList<PersonIndex> Index) {
		double[][] sc = new double[Index.size()][1];
		String key;
		PersonIndex pi;
		for (int i = 0; i < Index.size(); i++) {
			pi = Index.get(i);
			key = pi.getKey();
			sc[i][0] = ((Double) ScoreHashTable.get(key)).doubleValue();
		}
		return sc;
	}

	public void RabinowitzApproach() {
		if (PedFileType) {
			PedData.RabinowitzApproach();
		}
	}

	/**
	 * It is used for real data analysis.
	 */
	public void realCreateTable() {
		clearTables();
		if (Matched == null || unMatched == null) {
			return;
		}
		Enumeration famstrList = PedData.getFamStrList();
		String[] FID = new String[PedData.getNumFamilies()];
		int ind = 0;
		while (famstrList.hasMoreElements()) {
			FID[ind++] = (String) famstrList.nextElement();
		}
		Arrays.sort(FID);

		for (int i = 0; i < FID.length; i++) {
			String fid = FID[i];
			FamilyUnit FamUnit = PhenoData.getFamilyUnit(fid);
			FamilyStruct FamStr = PedData.getFamilyStruct(fid);
			Enumeration sibList = FamStr.getPersonList();
			String[] PID = new String[FamStr.getNumPersons()];
			int c = 0;
			while (sibList.hasMoreElements()) {
				PID[c++] = (String) sibList.nextElement();
			}
			Arrays.sort(PID);

			String pid;
			Person per;
			PseudoPerson pseudoper;
			Subject sub;
			for (int j = 0; j < PID.length; j++) {
				pid = PID[j];
				if (IsUnMatchedIndividual(fid, pid)) {
					continue;
				}
				try {
					per = FamStr.getPerson(pid);
					pseudoper = FamStr.getPseudoPerson(pid);
					sub = FamUnit.getSubject(pid);
					FullPersonTable.add(new PersonIndex(fid, pid));
					FullCovariateTable.add((ArrayList) sub.getTraits().clone());
					FullStatusTable.add(new Integer(per.getAffectedStatus()));
					if (!FamStr.hasAncestor(per.getPersonID())) {
						continue;
					}
					TDTGenoTable.add(per.getGenotype().clone());
					TransmittedTable.add(per.getGenotype().clone());
					if (PedFileType) {
						TDTGenoTable.add(pseudoper.getPseudoGenotype().clone());
						NontransmittedTable.add(pseudoper.getPseudoGenotype()
								.clone());
					}
					StatusTable.add(new Integer(per.getAffectedStatus()));
					CovariateTable.add(sub.getTraits().clone());
					PersonTable.add(new PersonIndex(fid, pid));
				} catch (MDRPedFileException e) {
					System.err.println("Can't find the individual " + pid
							+ " genotype in family " + fid);
				} catch (GMDRPhenoFileException e) {
					System.err.println("Can't find the individual " + pid
							+ " phenotype in family " + fid);
				}
			}
		}
	}

	/**
	 * It is used for real data analysis
	 * 
	 * @param TDTped
	 *            The file for printing out transmitted-nontransmitted
	 *            genotypes.
	 * @param TDTphe
	 *            The file for printing out score.
	 * @param TDTpi
	 *            The file for printing index of the genotypes.
	 * @param datatype
	 *            Determines which part of data to print out
	 * @throws IOException
	 * @throws CalEngineException
	 */
	public void realPrintGMDR(String TDTped, String TDTphe, String TDTpi,
			int datatype) throws IOException, CalEngineException {
		PrintWriter TDTpedout = new PrintWriter(TDTped);
		PrintWriter TDTpheout = new PrintWriter(TDTphe);
		PrintWriter TDTpiout = new PrintWriter(TDTpi);
		Hashtable faminformative = PedData.getFamInformative();
		PersonIndex PI;
		int countmiss = 0;
		Hashtable InformativePersonHash = new Hashtable();

		DiscardedIndividual = new ArrayList();
		CardedIndividual = new ArrayList();
		for (int i = 0; i < PersonTable.size(); i++) {
			PI = PersonTable.get(i);
			if ((!((Boolean) faminformative.get(PI.getFamilyID()))
					.booleanValue())
					|| (!ScoreHashTable.containsKey(PI.getKey())))// not
			// informative
			// or no
			// score
			{
				countmiss++;
				DiscardedIndividual.add(PI.getKey());
				// System.out.println(PI.getKey());
				continue;
			} else {
				CardedIndividual.add(PI);
			}
			InformativePersonHash.put(new String(PI.getKey()), new Integer(i));
			TDTpiout.println(PI.getFamilyID() + "\t" + PI.getIndividualID());
		}
		TDTpiout.close();
		if (countmiss > 0) {
			System.out.println("=================");
			System.out.println("Discarded " + countmiss + " individuals.");
			for (int i = 0; i < DiscardedIndividual.size(); i++) {
				System.out.println((String) DiscardedIndividual.get(i));
			}
			System.out.println("=================");
		}
		ArrayList marker;

		ArrayList markerinfor = PedData.getMarkerInformation();
		for (int i = 0; i < markerinfor.size(); i++) {
			TDTpedout.print(markerinfor.get(i) + "\t");
		}
		TDTpedout.println("class");

		int sindex = 0;
		ArrayList<String> PrintOutPersonTable = new ArrayList();
		for (int i = 0; i < TransmittedTable.size(); i++) {
			PI = PersonTable.get(i);
			if (!InformativePersonHash.containsKey(PI.getKey())) {
				continue;
			}
			if ((((Integer) StatusTable.get(i)).intValue() - 1) == 1
					&& datatype == 1) {
				continue;
			}

			marker = (ArrayList) TransmittedTable.get(i);
			for (int j = 0; j < marker.size(); j++) {
				TDTpedout.print(marker.get(j) + "\t");
			}
			TDTpedout.println(((Integer) StatusTable.get(i)).intValue() - 1);
			marker = (ArrayList) NontransmittedTable.get(i);
			for (int j = 0; j < marker.size(); j++) {
				TDTpedout.print(marker.get(j) + "\t");
			}
			TDTpedout
					.println(1 - (((Integer) StatusTable.get(i)).intValue() - 1));
			sindex++;
			PrintOutPersonTable.add(new String(PI.getKey()));
		}
		TDTpedout.close();
		try {
			if (PheIndex[0] == -1) {
				TDTpheout.println("Affection_status");
			} else {
				TDTpheout.println(PhenoData.getTraitAtI(PheIndex[0]));
			}
		} catch (GMDRPhenoFileException e) {
			System.err.println("Could not find the phenotype at index "
					+ PheIndex[0]);
		}
		for (int i = 0; i < PrintOutPersonTable.size(); i++) {
			String pi = PrintOutPersonTable.get(i);
			double s = ((Double) ScoreHashTable.get(pi)).doubleValue();
			TDTpheout.println(s);
			TDTpheout.println((-1) * s);
		}
		TDTpheout.close();
	}

	public void PrintNullGMDR(String TDTped, int datatype, long seed)
			throws IOException, CalEngineException {
		PrintWriter TDTpedout = new PrintWriter(TDTped);
		Hashtable faminformative = PedData.getFamInformative();
		PersonIndex PI;
		int countmiss = 0;
		Hashtable InformativePersonHash = new Hashtable();

		for (int i = 0; i < PersonTable.size(); i++) {
			PI = PersonTable.get(i);
			if ((!((Boolean) faminformative.get(PI.getFamilyID()))
					.booleanValue())
					|| (!ScoreHashTable.containsKey(PI.getKey()))){
				continue;
			}
			InformativePersonHash.put(new String(PI.getKey()), new Integer(i));
		}

		ArrayList marker;
		ArrayList markerinfor = PedData.getMarkerInformation();
		for (int i = 0; i < markerinfor.size(); i++) {
			TDTpedout.print(markerinfor.get(i) + "\t");
		}
		TDTpedout.println("class");

		Random rnd = new Random(seed);
		int PI_fam_idx = -1;
		double r = 0;
		for (int i = 0; i < TransmittedTable.size(); i++) {
			PI = PersonTable.get(i);
			if (!InformativePersonHash.containsKey(PI.getKey())) {
				continue;
			}
			if ((((Integer) StatusTable.get(i)).intValue() - 1) == 1
					&& datatype == 1) {
				continue;
			}
			if((Integer.parseInt(PI.FamilyID)) != PI_fam_idx) {//switch for a new family
				r = rnd.nextDouble();
				PI_fam_idx = Integer.parseInt(PI.FamilyID);
			}
			if (r > 0.5) {
				marker = (ArrayList) TransmittedTable.get(i);
				for (int j = 0; j < marker.size(); j++) {
					TDTpedout.print(marker.get(j) + "\t");
				}
				TDTpedout.println(((Integer) StatusTable.get(i)).intValue() - 1);
				marker = (ArrayList) NontransmittedTable.get(i);
				for (int j = 0; j < marker.size(); j++) {
					TDTpedout.print(marker.get(j) + "\t");
				}
				TDTpedout
					.println(1 - (((Integer) StatusTable.get(i)).intValue() - 1));
			} else {//switch genotypes
				marker = (ArrayList) NontransmittedTable.get(i);
				for (int j = 0; j < marker.size(); j++) {
					TDTpedout.print(marker.get(j) + "\t");
				}
				TDTpedout.println(((Integer) StatusTable.get(i)).intValue() - 1);
				marker = (ArrayList) TransmittedTable.get(i);
				for (int j = 0; j < marker.size(); j++) {
					TDTpedout.print(marker.get(j) + "\t");
				}
				TDTpedout
					.println(1 - (((Integer) StatusTable.get(i)).intValue() - 1));
			}
		}
		TDTpedout.close();
	}

	/**
	 * It is used for real data analysis
	 * 
	 * @param TDTped
	 *            The file for printing out transmitted-nontransmitted
	 *            genotypes.
	 * @param TDTphe
	 *            The file for printing out score.
	 * @param TDTpi
	 *            The file for printing index of the genotypes.
	 * @param datatype
	 *            Determines which part of data to print out
	 * @throws IOException
	 * @throws CalEngineException
	 */
	public void RabinowitzPrintGMDR(String TDTped, String TDTphe,
			boolean isPermutation) throws IOException, CalEngineException {
		// datatype = 1 for transmitted only or else untransmitted.
		PrintWriter TDTpedout = new PrintWriter(TDTped);
		PrintWriter TDTpheout = null;
		if (!isPermutation) {
			TDTpheout = new PrintWriter(TDTphe);
		}
		Hashtable faminformative = PedData.getFamInformative();
		PersonIndex PI;
		int countmiss = 0;
		Hashtable InformativePersonHash = new Hashtable();

		ArrayList Discard = new ArrayList();
		for (int i = 0; i < PersonTable.size(); i++) {
			PI = PersonTable.get(i);
			if ((!((Boolean) faminformative.get(PI.getFamilyID()))
					.booleanValue())
					|| (!ScoreHashTable.containsKey(PI.getKey()))) { // not
				// informative
				// or no
				// score
				countmiss++;
				Discard.add(PI.getKey());
				continue;
			}
			InformativePersonHash.put(new String(PI.getKey()), new Integer(i));
		}
		if (countmiss > 0) {
			System.out.println("=================");
			System.out.println("Discarded " + countmiss + " individuals.");
			for (int i = 0; i < Discard.size(); i++) {
				System.out.println((String) Discard.get(i));
			}
			System.out.println("=================");
		}
		ArrayList marker;

		ArrayList markerinfor = PedData.getMarkerInformation();
		for (int i = 0; i < markerinfor.size(); i++) {
			TDTpedout.print(markerinfor.get(i) + "\t");
		}
		TDTpedout.println("class");

		ArrayList<String> PrintOutPersonTable = new ArrayList();
		for (int i = 0; i < TransmittedTable.size(); i++) {
			PI = PersonTable.get(i);
			if (!InformativePersonHash.containsKey(PI.getKey())) {
				continue;
			}
			PrintOutPersonTable.add(new String(PI.getKey()));
			marker = (ArrayList) TransmittedTable.get(i);
			if (!isPermutation) {
				for (int j = 0; j < marker.size(); j++) {
					TDTpedout.print(marker.get(j) + "\t");
				}
				TDTpedout
						.println(((Integer) StatusTable.get(i)).intValue() - 1);
			} else {
				marker = (ArrayList) RabinowitzTable.get(i);
				for (int j = 0; j < marker.size(); j++) {
					TDTpedout.print(marker.get(j) + "\t");
				}
				TDTpedout
						.println(((Integer) StatusTable.get(i)).intValue() - 1);
			}
		}
		TDTpedout.close();
		if (!isPermutation) {
			try {
				TDTpheout.println(PhenoData.getTraitAtI(PheIndex[0]));
			} catch (GMDRPhenoFileException e) {
				System.err.println("Could not find the phenotype at index "
						+ PheIndex[0]);
			}
			for (int i = 0; i < PrintOutPersonTable.size(); i++) {
				String pi = PrintOutPersonTable.get(i);
				double s = ((Double) ScoreHashTable.get(pi)).doubleValue();
				TDTpheout.println(s);// print the scores for transtmitted
				// phenotype only;
			}
			TDTpheout.close();
		}
	}

	private void setAdjustmentOfCovariate(boolean adj) {
		AdjustmentOfCovariate = adj;
	}

	public void setFileType(boolean filetype) {
		PedFileType = filetype;
	}

	private void setSelectedCovariateIndex(int[] CIndex) {
		if (CIndex != null) {
			CovIndex = new int[CIndex.length];
			System.arraycopy(CIndex, 0, CovIndex, 0, CIndex.length);
		} else {
			CovIndex = null;
		}
	}

	private void setScoreProc(int sp) {
		scoreProc = sp;
	}

	private void setSelectedPhenoIndex(int[] PIndex) {
		PheIndex = new int[PIndex.length];
		System.arraycopy(PIndex, 0, PheIndex, 0, PIndex.length);
	}

	public void setTDTScore(double[][] score) {
		TDTScore = new double[score.length][score[0].length];
		for (int i = 0; i < score.length; i++) {
			for (int j = 0; j < score[i].length; j++) {
				TDTScore[i][j] = score[i][j];
			}
		}
	}

	private void TableConversion() {
		Status = new double[StatusTable.size()][1];
		Covariates = new double[CovariateTable.size()][PhenoData.getNumTraits()];
		FullStatus = new double[FullStatusTable.size()][1];
		FullCovariates = new double[FullCovariateTable.size()][PhenoData
				.getNumTraits()];
		ArrayList cov;
		String co;
		for (int i = 0; i < StatusTable.size(); i++) {
			Status[i][0] = ((Integer) StatusTable.get(i)).intValue();
			cov = (ArrayList) CovariateTable.get(i);
			for (int j = 0; j < cov.size(); j++) {
				co = (String) cov.get(j);
				Covariates[i][j] = Double.parseDouble(co);
			}
		}

		for (int i = 0; i < FullStatusTable.size(); i++) {
			FullStatus[i][0] = ((Integer) FullStatusTable.get(i)).intValue();
			cov = (ArrayList) FullCovariateTable.get(i);
			for (int j = 0; j < cov.size(); j++) {
				co = (String) cov.get(j);
				FullCovariates[i][j] = Double.parseDouble(co);
			}
		}
	}
}
