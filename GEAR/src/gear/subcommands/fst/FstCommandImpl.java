package gear.subcommands.fst;

import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import gear.data.InputDataSet2;
import gear.data.SubjectID;
import gear.family.GenoMatrix.GenotypeMatrix;
import gear.family.pedigree.file.SNP;
import gear.family.plink.PLINKParser;
import gear.family.qc.rowqc.SampleFilter;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.BufferedReader;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;

public class FstCommandImpl extends CommandImpl {

	private FstCommandArguments fstArgs;
	private SampleFilter sf;
	private GenotypeMatrix pGM;

	private InputDataSet2 data = new InputDataSet2();
	private HashMap<SubjectID, Integer> fstGroup = NewIt.newHashMap();
	private int fstG = 0;
	private double[][] Fst = null;

	@Override
	public void execute(CommandArguments cmdArgs) {
		fstArgs = (FstCommandArguments) cmdArgs;
		PLINKParser pp = PLINKParser.parse(cmdArgs);

		data.addFile(fstArgs.getFam()); // geno
		data.addFile(fstArgs.getGroupFile()); // group
		if (fstArgs.isKeepFile())
			data.addFile(fstArgs.getKeepFile()); // keep
		data.LineUpFiles();

		sf = new SampleFilter(pp.getPedigreeData(), data.getMatchSubjetList());
		pGM = new GenotypeMatrix(sf.getSample(), pp.getMapData(), cmdArgs);

		readGroup();
		calFst();
		writeFst();
	}

	private void readGroup() {

		BufferedReader reader = BufferedReader.openTextFile(fstArgs.getGroupFile(), "Fst::group file");
		String[] tokens = reader.readTokensAtLeast(3);

		if (tokens == null) {
			Logger.printUserError("The group file '" + fstArgs.getGroupFile() + "' is empty.");
			System.exit(1);
		}

		int numCols = tokens.length;
		HashMap<String, Integer> groupCntMap = NewIt.newHashMap();
		ArrayList<String> labelList = NewIt.newArrayList();

		do {
			SubjectID subjectID = new SubjectID(/* famID */tokens[0], /* indID */tokens[1]);
			if (fstGroup.containsKey(subjectID)) {
				reader.errorPreviousLine("Subject " + subjectID + " is duplicated.");
			}
			if (!groupCntMap.containsKey(tokens[2])) {
				groupCntMap.put(tokens[2], 1);
				fstGroup.put(subjectID, 1);
				labelList.add(tokens[2]);
			} else {
				for (Iterator<String> i = groupCntMap.keySet().iterator(); i.hasNext();) {
					String label = i.next();
					if (label.equals((String) tokens[2])) {
						int cn = groupCntMap.get(label).intValue();
						cn++;
						groupCntMap.put(label, cn);
					}
				}
			}

			int cnt = 0;
			for (Iterator<String> i = labelList.iterator(); i.hasNext();) {
				String label = i.next();
				if (label.equals((String) tokens[2])) {
					fstGroup.put(subjectID, cnt);
				}
				cnt++;
			}

		} while ((tokens = reader.readTokens(numCols)) != null);
		reader.close();

		if (groupCntMap.size() < 2) {
			Logger.printUserLog("Should be at least 2 groups.");
			System.exit(1);
		}
		Logger.printUserLog("Read " + fstGroup.size() + " subjects for " + groupCntMap.size() + " groups.");
		fstG = groupCntMap.size();

		for (Iterator<String> i = groupCntMap.keySet().iterator(); i.hasNext();) {
			String label = i.next();
			Logger.printUserLog(groupCntMap.get(label) + " subjects are grouped to '" + label + "'.");
		}
	}

	private void calFst() {

		Fst = new double[pGM.getNumMarker()][2];

		for (int i = 0; i < pGM.getNumMarker(); i++) {
			double[] gSub = new double[fstG], freqSub = new double[fstG];
			double[] nSub = new double[fstG];
			double[] hetSub = new double[fstG], hetFreq = new double[fstG];

			for (int j = 0; j < pGM.getNumIndivdial(); j++) {
				SubjectID sub = new SubjectID(pGM.getSample().get(j).getFamilyID(),
						pGM.getSample().get(j).getIndividualID());
				int gIdx = fstGroup.get(sub).intValue();
				if (!pGM.isNA(j, i)) {
					nSub[gIdx] += 1.0;
					gSub[gIdx] += pGM.getAdditiveScore(j, i) * 1.0;
					if (pGM.getAdditiveScore(j, i) == 1) {
						hetSub[gIdx] += 1.0;
					}
				}
			}

			boolean badLocus = false;
			double nTotal = 0, nMean = 0, hetTotal = 0;
			double gTotal = 0, freqMean = 0, hetMean = 0;

			for (int k = 0; k < fstG; k++) {
				if (nSub[k] == 0) {
					badLocus = true;
					break;
				}
				nTotal += nSub[k];
				freqSub[k] = gSub[k] / (2.0 * nSub[k]);
				hetFreq[k] = hetSub[k] / (nSub[k]);
				hetTotal += hetSub[k];
				gTotal += gSub[k];
			}
			// see bruce weir [genetic data analysis] page 172-3.
			if (!badLocus) {
				nMean = nTotal / (1.0 * fstG);
				freqMean = gTotal / (2.0 * nTotal);
				hetMean = hetTotal / (1.0 * nTotal);
				double nc = nTotal / (fstG - 1);
				double SAsq = 0;
				for (int k = 0; k < fstG; k++) {
					nc -= (nSub[k] * nSub[k]) / nTotal;
					SAsq += nSub[k] * (freqSub[k] - freqMean) * (freqSub[k] - freqMean) / ((fstG - 1) * nMean);
				}
				nc /= (fstG - 1);

				double T1 = 0, T2 = 0;
				T1 = SAsq - 1 / (nMean - 1) * (freqMean * (1 - freqMean) - SAsq * (fstG - 1) / fstG);
				T2 = (nc - 1) / (nMean - 1) * freqMean * (1 - freqMean)
						+ (1 + (fstG - 1) * (nMean - nc) / (nMean - 1)) * (SAsq / fstG);
				Fst[i][0] = T1 / T2;
				Fst[i][1] = SAsq / (freqMean * (1 - freqMean));
				// SNP snp = pGM.getSNPList().get(i);
				// Logger.printUserLog("SNP::" + snp.getName() + " FstR " + Fst[i][0] + " FstF "
				// + Fst[i][1]+ " Freq " + freqMean + " " + SAsq + " "
				// + T1 + " " + T2);
			}
		}
	}

	private void writeFst() {
		PrintStream fstOut = FileUtil.CreatePrintStream(this.fstArgs.getOutRoot() + ".fst");
		DecimalFormat df = new DecimalFormat("0.0000");
		DecimalFormat dfE = new DecimalFormat("0.00E000");

		double s1 = 0, s2 = 0;
		double ss1 = 0, ss2 = 0;
		int cnt1 = 0, cnt2 = 0;
		fstOut.println("SNP\tCHR\tBP\tRefAllele\tAltAllele\tFst[R]\tFst[F]");
		for (int i = 0; i < pGM.getNumMarker(); i++) {
			SNP snp = pGM.getSNPList().get(i);
			fstOut.print(snp.getName() + "\t" + snp.getChromosome() + "\t" + snp.getPosition() + "\t"
					+ snp.getFirstAllele() + "\t" + snp.getSecAllele() + "\t");
			if (Math.abs(Fst[i][0]) > 0.0001) {
				fstOut.print(df.format(Fst[i][0]) + "\t");
			} else {
				fstOut.print(dfE.format(Fst[i][0]) + "\t");
			}
			if (Math.abs(Fst[i][1]) < 0.0001) {
				fstOut.println(df.format(Fst[i][1]));
			} else {
				fstOut.println(dfE.format(Fst[i][1]));
			}

			if (Fst[i][0] != 0) {
				ss1 += Fst[i][0] * Fst[i][0];
				s1 += Fst[i][0];
				cnt1++;
			}
			if (Fst[i][1] != 0) {
				ss2 += Fst[i][1] * Fst[i][1];
				s2 += Fst[i][1];
				cnt2++;
			}
		}
		fstOut.close();
		Logger.printUserLog(
				"Fst[R] mean: " + s1 / (1.0 * cnt1) + ", sd: " + Math.sqrt((ss1/cnt1 - s1 * s1 / (cnt1 * cnt1))));
		Logger.printUserLog(
				"Fst[F] mean: " + s2 / (1.0 * cnt2) + ", sd: " + Math.sqrt((ss2/cnt1 - s2 * s2 / (cnt2 * cnt2))));
	}
}
