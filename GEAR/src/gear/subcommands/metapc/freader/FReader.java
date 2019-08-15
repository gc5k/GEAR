package gear.subcommands.metapc.freader;

import gear.ConstValues;
import gear.util.BufferedReader;
import gear.util.Logger;
import gear.util.NewIt;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

public class FReader {

	public FReader(String[] MetaFile, boolean[] FileKeep, String[] field, boolean isQT, boolean isGZ, boolean isChr,
			HashSet<String> Chr) {
		if(field != null ) {
			this.field = field;
		}
		this.isGZ = isGZ;

		workingMetaFile = Collections.synchronizedList(NewIt.newArrayList());
		this.isChrKeep = isChr;
		this.chrSet = Chr;

		for (int i = 0; i < FileKeep.length; i++) {
			if (FileKeep[i]) {
				workingMetaFile.add(MetaFile[i]);
			}
		}

		if (workingMetaFile.size() == 0) {
			Logger.printUserLog("No cohorts left for analysis. GEAR quit.");
			System.exit(0);
		} else {
			Logger.printUserLog(workingMetaFile.size() + " cohorts are remained for analysis.");
		}

		KeyIdx = new int[workingMetaFile.size()][FConstant.num_key];
		for (int i = 0; i < KeyIdx.length; i++) {
			Arrays.fill(KeyIdx[i], -1);
		}
	}

	public void Start() {
		for (int i = 0; i < workingMetaFile.size(); i++) {
			HashMap<String, FStat> m = readMeta(i);
			MStat.add(m);
		}
	}

	public int getCohortNum() {
		return workingMetaFile.size();
	}

	public String[] getMetaFile() {
		return workingMetaFile.toArray(new String[0]);
	}

	public int[][] getKeyIndex() {
		return KeyIdx;
	}

	public int getNumMetaFile() {
		return workingMetaFile.size();
	}

	public ArrayList<HashMap<String, FStat>> getMetaStat() {
		return MStat;
	}

	public List<ArrayList<String>> getMetaSNPArray() {
		return MetaSNPArray;
	}

	public HashMap<String, ArrayList<Integer>> getMetaSNPTable() {
		return MetaSNPTable;
	}

	private HashMap<String, FStat> readMeta(int metaIdx) {
		BufferedReader reader = null;
		if (isGZ) {
			reader = BufferedReader.openGZipFile(workingMetaFile.get(metaIdx), "Summary Statistic file");
		} else {
			reader = BufferedReader.openTextFile(workingMetaFile.get(metaIdx), "Summary Statistic file");
		}
		Logger.printUserLog("");
		Logger.printUserLog("Reading summary stats from '" + workingMetaFile.get(metaIdx) + "'...");

		String[] tokens = reader.readTokens();
		int tokenLen = tokens.length;

		for (int i = 0; i < tokens.length; i++) {
			if (tokens[i].equalsIgnoreCase(field[FConstant.SNP])) {
				KeyIdx[metaIdx][FConstant.SNP] = i;
			}
			if (tokens[i].equalsIgnoreCase(field[FConstant.CHR])) {
				KeyIdx[metaIdx][FConstant.CHR] = i;
			}
			if (tokens[i].equalsIgnoreCase(field[FConstant.A1])) {
				KeyIdx[metaIdx][FConstant.A1] = i;
			}
			if (tokens[i].equalsIgnoreCase(field[FConstant.A2])) {
				KeyIdx[metaIdx][FConstant.A2] = i;
			}
			if (tokens[i].equalsIgnoreCase(field[FConstant.Fvalue])) {
				KeyIdx[metaIdx][FConstant.Fvalue] = i;
			}
			// if (tokens[i].equalsIgnoreCase(field[FConstant.BP]))
			// {
			// KeyIdx[metaIdx][FConstant.BP] = i;
			// }
			// if (tokens[i].equalsIgnoreCase(field[FConstant.N]))
			// {
			// KeyIdx[metaIdx][FConstant.N] = i;
			// }
		}

		boolean qFlag = false;

		if (KeyIdx[metaIdx][FConstant.SNP] == -1) {
			Logger.printUserLog("Cannot find the snp column in " + workingMetaFile.get(metaIdx));
			qFlag = true;
		}

		if (KeyIdx[metaIdx][FConstant.CHR] == -1) {
			Logger.printUserLog("Cannot find the chr column in " + workingMetaFile.get(metaIdx));
			qFlag = true;
		}
		/*
		 * if (KeyIdx[metaIdx][FConstant.BP] == -1) {
		 * Logger.printUserLog("Cannot find the bp column in " +
		 * workingMetaFile.get(metaIdx)); }
		 */
		if (KeyIdx[metaIdx][FConstant.Fvalue] == -1) {
			Logger.printUserLog("Cannot find the value column in " + workingMetaFile.get(metaIdx));
			qFlag = true;
		}

		if (KeyIdx[metaIdx][FConstant.A1] == -1) {
			Logger.printUserLog("Cannot find the allele 1 column in " + workingMetaFile.get(metaIdx));
			qFlag = true;
		}

		if (KeyIdx[metaIdx][FConstant.A2] == -1) {
			Logger.printUserLog("Cannot find the allele 2 column in " + workingMetaFile.get(metaIdx));
			qFlag = true;
		}

		if (qFlag) {
			Logger.printUserLog("GEAR quit.");
			System.exit(0);
		}

		HashMap<String, FStat> sumstat = NewIt.newHashMap();
		ArrayList<String> snpArray = NewIt.newArrayList();
		int total = 0;
		int cnt = 0;
		int cntDup = 0;
		int cntBadChr = 0;
		// int cntBadBp = 0;
		int cntBadFvalue = 0;
		int cntBadA1 = 0;
		int cntBadA2 = 0;
		int cntMissSNP = 0;

		while ((tokens = reader.readTokens(tokenLen)) != null) {
			total++;
			if (FConstant.isNASNP(tokens[KeyIdx[metaIdx][FConstant.SNP]])) {
				cntMissSNP++;
				continue;
			}

			if (ConstValues.isNA(tokens[KeyIdx[metaIdx][FConstant.CHR]])) {
				cntBadChr++;
				continue;
			}

			// if (KeyIdx[metaIdx][FConstant.BP] != -1 &&
			// ConstValues.isNA(tokens[KeyIdx[metaIdx][FConstant.BP]]))
			// {
			// cntBadBp++;
			// continue;
			// }

			if (ConstValues.isNA(tokens[KeyIdx[metaIdx][FConstant.Fvalue]])) {
				cntBadFvalue++;
				continue;
			}

			if (tokens[KeyIdx[metaIdx][FConstant.A1]].length() != 1) {
				cntBadA1++;
				continue;
			}

			if (KeyIdx[metaIdx][FConstant.A2] != -1 && tokens[KeyIdx[metaIdx][FConstant.A2]].length() != 1) {
				cntBadA2++;
				continue;
			}

			FStat ms = new FStat(tokens[KeyIdx[metaIdx][FConstant.SNP]],
					Float.parseFloat(tokens[KeyIdx[metaIdx][FConstant.Fvalue]]),
					tokens[KeyIdx[metaIdx][FConstant.A1]].charAt(0), tokens[KeyIdx[metaIdx][FConstant.A2]].charAt(0));
			if (KeyIdx[metaIdx][FConstant.CHR] != -1) {
				if (tokens[KeyIdx[metaIdx][FConstant.CHR]].equalsIgnoreCase("X")) {
					ms.setChr(23);
				} else if (tokens[KeyIdx[metaIdx][FConstant.CHR]].equalsIgnoreCase("Y")) {
					ms.setChr(24);
				} else if (tokens[KeyIdx[metaIdx][FConstant.CHR]].equalsIgnoreCase("XY")) {
					ms.setChr(25);
				} else if (tokens[KeyIdx[metaIdx][FConstant.CHR]].equalsIgnoreCase("MT")) {
					ms.setChr(26);
				} else {
					int chr = -1;
					try {
						chr = Integer.parseInt(tokens[KeyIdx[metaIdx][FConstant.CHR]]);
					} catch (NumberFormatException e) {
						Logger.printUserLog(e.toString() + " in line " + total + " in '" + workingMetaFile.get(metaIdx)
								+ ".' is a bad value for chromosome. Skipped this marker.");
						continue;
					}
					ms.setChr(chr);
				}

				if (isChrKeep) {
					if (chrSet != null && !chrSet.contains((new Integer(ms.getChr())).toString())) {
						continue;
					}
				} else {
					if (chrSet != null && chrSet.contains((new Integer(ms.getChr())).toString())) {
						continue;
					}
				}
			}
			// if (KeyIdx[metaIdx][FConstant.BP] != -1)
			// {
			// long bp = -1;
			// try
			// {
			// bp = Long.parseLong(tokens[KeyIdx[metaIdx][FConstant.BP]]);
			// }
			// catch (NumberFormatException e)
			// {
			// Logger.printUserLog(e.toString() + " in line " + total + " in '" +
			// workingMetaFile.get(metaIdx) + "' is a bad value for position. Skipped this
			// marker.");
			// continue;
			// }
			// ms.setBP(bp);
			// }

			if (sumstat.containsKey(ms.getSNP())) {
				Logger.printUserLog(
						"Warning: Marker '" + ms.getSNP() + "' duplicated input, first instance used, others skipped.");
				cntDup++;
			} else {
				sumstat.put(ms.getSNP(), ms);
				snpArray.add(ms.getSNP());

				if (MetaSNPTable.containsKey(ms.getSNP())) {
					ArrayList<Integer> snpCnt = MetaSNPTable.get(ms.getSNP());
					snpCnt.set(metaIdx, 1);
					Integer Int = snpCnt.get(snpCnt.size() - 1);
					Int++;
					snpCnt.set(snpCnt.size() - 1, Int);
				} else {
					ArrayList<Integer> snpCnt = NewIt.newArrayList();
					snpCnt.ensureCapacity(workingMetaFile.size() + 1);
					for (int ii = 0; ii < workingMetaFile.size() + 1; ii++) {
						snpCnt.add(0);
					}
					snpCnt.set(metaIdx, 1);
					snpCnt.set(snpCnt.size() - 1, 1);
					MetaSNPTable.put(ms.getSNP(), snpCnt);
				}
				cnt++;
			}
		}
		reader.close();

		if (cntMissSNP > 0) {
			String lc = cntMissSNP == 1 ? "locus" : "loci";
			Logger.printUserLog("Removed " + cntMissSNP + " " + lc + " due to bad marker name(s)");
		}

		if (cntBadChr > 0) {
			String lc = cntBadChr == 1 ? "locus" : "loci";
			Logger.printUserLog("Removed " + cntBadChr + " " + lc + " due to incorrect chromosome(s).");
		}

		// if (cntBadBp > 0)
		// {
		// String lc = cntBadBp == 1 ? "locus" : "loci";
		// Logger.printUserLog("Removed " + cntBadBp + " " + lc + " due to incorrect
		// physical position(s).");
		// }

		if (cntBadFvalue > 0) {
			String lc = cntBadFvalue == 1 ? "locus" : "loci";
			Logger.printUserLog("Removed " + cntBadFvalue + " " + lc + " due to incorrect effect(s).");
		}

		if (cntBadFvalue > 0) {
			String lc = cntBadFvalue == 1 ? "locus" : "loci";
			Logger.printUserLog("Removed " + cntBadFvalue + " " + lc + " due to incorrect se.");
		}

		if (cntBadA1 > 0) {
			String lc = cntBadA1 == 1 ? "locus" : "loci";
			Logger.printUserLog("Removed " + cntBadA1 + " " + lc + " due to bad a1 allele(s).");
		}

		if (cntBadA2 > 0) {
			String lc = cntBadA2 == 1 ? "locus" : "loci";
			Logger.printUserLog("Removed " + cntBadA2 + " " + lc + " due to bad a2 allele(s).");
		}

		if (cntDup > 0) {
			String lc = cntDup == 1 ? "locus" : "loci";
			Logger.printUserLog("Removed " + cntDup + " duplicated " + lc);
		}

		MetaSNPArray.add(snpArray);
		if (cnt == 0) {
			Logger.printUserLog("Did not find any summary statistics from '" + workingMetaFile.get(metaIdx) + ".'");
			System.exit(0);
		} else {
			Logger.printUserLog("Read " + cnt + " (of " + total + ") summary statistics from '"
					+ workingMetaFile.get(metaIdx) + "'.");
		}

		return sumstat;
	}

	public List<String> getWorkingMetaFile() {
		return workingMetaFile;
	}

	private boolean isChrKeep = true;

	private HashSet<String> chrSet = null;

	private String[] field = {"SNP", "CHR", "A1", "A2", "MAF"};
	private boolean isGZ;
	private int[][] KeyIdx; // snp, chr, bp, beta, se, p, a1, a2
	private List<String> workingMetaFile;

	private ArrayList<HashMap<String, FStat>> MStat = NewIt.newArrayList();

	private List<ArrayList<String>> MetaSNPArray = Collections.synchronizedList(NewIt.newArrayList());
//	private ArrayList<ArrayList<String>> MetaSNPArray = NewIt.newArrayList();
	private HashMap<String, ArrayList<Integer>> MetaSNPTable = NewIt.newHashMap();

	// private int[] keepCohortIdx;

}
