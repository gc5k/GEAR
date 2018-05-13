package gear.subcommands.grm;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;

import gear.data.Person;
import gear.family.GenoMatrix.GenotypeMatrix;
import gear.family.pedigree.PersonIndex;
import gear.family.plink.PLINKParser;
import gear.family.qc.rowqc.SampleFilter;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.pop.PopStat;

public class GRMCommandImpl extends CommandImpl {
	private GenotypeMatrix pGM;

	private double[][] allelefreq;
	private double[] allelevar;
	private GRMCommandArguments grmArgs;
	private SampleFilter sf;

	@Override
	public void execute(CommandArguments cmdArgs) {
		grmArgs = (GRMCommandArguments) cmdArgs;

		PLINKParser pp = PLINKParser.parse(grmArgs);
		sf = new SampleFilter(pp.getPedigreeData(), cmdArgs);
		pGM = new GenotypeMatrix(sf.getSample(), pp.getMapData(), cmdArgs);

		allelefreq = PopStat.calAlleleFrequency(pGM);
		allelevar = PopStat.calGenoVariance(pGM);

		makeAddScore();

		if (grmArgs.isDom()) {
			makeDomScore();
		}
	}

	public void makeDomScore() {
		double grmDomMean = 0;
		double grmDomSq = 0;

		StringBuffer sb = new StringBuffer();
		sb.append(grmArgs.getOutRoot());
		PrintStream grm = null;
		BufferedWriter grmGZ = null;
		if (grmArgs.isGZ()) {
			sb.append(".dom.grm.gz");
			grmGZ = FileUtil.ZipFileWriter(sb.toString());
		} else {
			sb.append(".dom.grm.txt");
			grm = FileUtil.CreatePrintStream(sb.toString());
		}

		int cnt = 0;
		for (int i = 0; i < pGM.getGRow(); i++) {
			for (int j = 0; j <= i; j++) {
				double[] d = GRMDomScore(i, j);
				if (i != j) {
					grmDomMean += d[1];
					grmDomSq += d[1] * d[1];
					cnt++;
				}
				if (grmArgs.isGZ()) {
					try {
						grmGZ.append((i + 1) + "\t" + (j + 1) + "\t" + d[0] + "\t" + d[1] + "\n");
					} catch (IOException e) {
						Logger.handleException(e,
								"error in writing '" + sb.toString() + "' for " + (i + 1) + " " + (j + 1) + ".");
					}
				} else {
					grm.println((i + 1) + "\t" + (j + 1) + "\t" + d[0] + "\t" + d[1]);
				}
			}
		}

		if (grmArgs.isGZ()) {
			try {
				grmGZ.close();
			} catch (IOException e) {
				Logger.handleException(e, " error in closing '" + sb.toString() + "'.");
			}
		} else {
			grm.close();
		}
		Logger.printUserLog("Writing GRM dominance scores into '" + sb.toString() + "'.");
		StringBuffer sb_id = new StringBuffer();
		sb_id.append(grmArgs.getOutRoot());
		PrintStream grm_id = null;
		sb_id.append(".dom.grm.id");
		grm_id = FileUtil.CreatePrintStream(sb_id.toString());

		ArrayList<PersonIndex> PI = sf.getSample();
		for (int i = 0; i < PI.size(); i++) {
			grm_id.println(PI.get(i).getFamilyID() + "\t" + PI.get(i).getIndividualID());
		}
		grm_id.close();
		Logger.printUserLog("Writing " + PI.size() + " individuals' information into '" + sb_id.toString() + "'.");

		grmDomMean /= cnt;
		grmDomSq /= cnt;
		double Effective_DomSample = -1 / grmDomMean;
		double grmDomSD = (grmDomSq - grmDomMean * grmDomMean) * cnt / (cnt - 1);
		double Effeictive_DomMarker = 1 / grmDomSD;

		DecimalFormat df = new DecimalFormat("0.0000");
		DecimalFormat dfE = new DecimalFormat("0.00E0");
		if (Math.abs(grmDomMean) > 0.0001) {
			Logger.printUserLog("Mean of dominance genetic relatedness is: " + df.format(grmDomMean));
		} else {
			Logger.printUserLog("Mean of dominance genetic relatedness is: " + dfE.format(grmDomMean));
		}

		if (Math.abs(grmDomSD) > 0.0001) {
			Logger.printUserLog("Sampling variance of dominance genetic relatedness is: " + df.format(grmDomSD));
		} else {
			Logger.printUserLog("Sampling variance of dominance genetic relatedness is: " + dfE.format(grmDomSD));
		}

		if (Math.abs(Effeictive_DomMarker) > 0.0001) {
			Logger.printUserLog("Effective dominance sample size is: " + df.format(Effective_DomSample));
		} else {
			Logger.printUserLog("Effective dominance sample size is: " + dfE.format(Effective_DomSample));
		}

		if (Math.abs(Effeictive_DomMarker) > 0.0001) {
			Logger.printUserLog("Effective number of dominance genome segments is: " + df.format(Effeictive_DomMarker));
		} else {
			Logger.printUserLog(
					"Effective number of dominance genome segments is: " + dfE.format(Effeictive_DomMarker));
		}
	}

	public void makeAddScore() {
		double grmMean = 0;
		double grmSq = 0;

		StringBuffer sb = new StringBuffer();
		sb.append(grmArgs.getOutRoot());
		PrintStream grm = null;
		BufferedWriter grmGZ = null;
		if (grmArgs.isGZ()) {
			sb.append(".grm.gz");
			grmGZ = FileUtil.ZipFileWriter(sb.toString());
		} else {
			sb.append(".grm.txt");
			grm = FileUtil.CreatePrintStream(sb.toString());
		}

		int cnt = 0;
		for (int i = 0; i < pGM.getGRow(); i++) {
			for (int j = 0; j <= i; j++) {
				double[] s = GRMScore(i, j);
				if (i != j) {
					grmMean += s[1];
					grmSq += s[1] * s[1];
					cnt++;
				}
				if (grmArgs.isGZ()) {
					try {
						grmGZ.append((i + 1) + "\t" + (j + 1) + "\t" + s[0] + "\t" + s[1] + "\n");
					} catch (IOException e) {
						Logger.handleException(e,
								"error in writing '" + sb.toString() + "' for " + (i + 1) + " " + (j + 1) + ".");
					}
				} else {
					grm.println((i + 1) + "\t" + (j + 1) + "\t" + s[0] + "\t" + s[1]);
				}
			}
		}

		if (grmArgs.isGZ()) {
			try {
				grmGZ.close();
			} catch (IOException e) {
				Logger.handleException(e, " error in closing '" + sb.toString() + "'.");
			}
		} else {
			grm.close();
		}
		Logger.printUserLog("Writing GRM scores into '" + sb.toString() + "'.");
		StringBuffer sb_id = new StringBuffer();
		sb_id.append(grmArgs.getOutRoot());
		PrintStream grm_id = null;
		sb_id.append(".grm.id");
		grm_id = FileUtil.CreatePrintStream(sb_id.toString());

		ArrayList<PersonIndex> PI = pGM.getSample();
		for (int i = 0; i < PI.size(); i++) {
			grm_id.println(PI.get(i).getFamilyID() + "\t" + PI.get(i).getIndividualID());
		}
		grm_id.close();
		Logger.printUserLog("Writing " + PI.size() + " individuals' information into '" + sb_id.toString() + "'.");

		grmMean /= cnt;
		grmSq /= cnt;
		double Effective_sample = -1 / grmMean + 1;
		double grmSD = (grmSq - grmMean * grmMean) * cnt / (cnt - 1);
		double Effeictive_marker = 1 / grmSD;

		DecimalFormat df = new DecimalFormat("0.0000");
		DecimalFormat dfE = new DecimalFormat("0.00E0");
		if (Math.abs(grmMean) > 0.0001) {
			Logger.printUserLog("Mean of genetic relatedness is: " + df.format(grmMean));
		} else {
			Logger.printUserLog("Mean of genetic relatedness is: " + dfE.format(grmMean));
		}

		if (Math.abs(grmSD) > 0.0001) {
			Logger.printUserLog("Sampling variance of genetic relatedness is: " + df.format(grmSD));
		} else {
			Logger.printUserLog("Sampling variance of genetic relatedness is: " + dfE.format(grmSD));
		}

		if (Math.abs(Effeictive_marker) > 0.0001) {
			Logger.printUserLog("Effective sample size is: " + df.format(Effective_sample));
		} else {
			Logger.printUserLog("Effective sample size is: " + dfE.format(Effective_sample));
		}

		if (Math.abs(Effeictive_marker) > 0.0001) {
			Logger.printUserLog("Effective number of genome segments is: " + df.format(Effeictive_marker));
		} else {
			Logger.printUserLog("Effective number of genome segments is: " + dfE.format(Effeictive_marker));
		}
	}

	private double[] GRMScore(int idx1, int idx2) {
		double[] s = { 0, 0 };

		for (int i = 0; i < allelefreq.length; i++) {
			
			//inside control
			if (allelefreq[i][1] == Double.NaN || allelefreq[i][0] == 0 || allelefreq[i][1] == 0 || allelevar[i] == 0) continue;

			int g1 = pGM.getAdditiveScore(idx1, i);
			int g2 = pGM.getAdditiveScore(idx2, i);
			double m = allelefreq[i][1];
			if (g1 == Person.MissingGenotypeCode || g2 == Person.MissingGenotypeCode) {
				continue;
			} else {
				double de = grmArgs.isInbred()? 4 * allelefreq[i][0] * allelefreq[i][1] : 2 * allelefreq[i][0] * allelefreq[i][1];
				if (grmArgs.isAdjVar()) de = allelevar[i];
				s[0]++;
				s[1] += (g1 - 2 * m) * (g2 - 2 * m) / de;
			}
		}

		if (s[0] > 0) {
			s[1] /= s[0];
		} else {
			s[0] = 0;
			s[1] = 0;
		}
		return s;
	}

	private double[] GRMDomScore(int idx1, int idx2) {
		double[] s = { 0, 0 };

		for (int i = 0; i < allelefreq.length; i++) {
			//inside control
			if (allelefreq[i][1] == Double.NaN || allelefreq[i][0] == 0 || allelefreq[i][1] == 0 || allelevar[i] == 0) continue;

			double m = allelefreq[i][1];
			int g1 = pGM.getAdditiveScore(idx1, i);
			int g2 = pGM.getAdditiveScore(idx2, i);
			if (g1 == Person.MissingGenotypeCode || g2 == Person.MissingGenotypeCode) {
				continue;
			} else {
				s[0]++;

				double s1, s2;
				if (g1 == 0) {
					s1 = 0;
				} else if (g1 == 1) {
					s1 = 2 * m;
				} else {
					s1 = (4 * m - 2);
				}

				if (g2 == 0) {
					s2 = 0;
				} else if (g2 == 1) {
					s2 = 2 * m;
				} else {
					s2 = (4 * m - 2);
				}

				s[1] += (s1 - 2 * m * m) * (s2 - 2 * m * m) / (4 * m * m * (1 - m) * (1 - m));
			}
		}

		if (s[0] > 0) {
			s[1] /= s[0];
		} else {
			s[0] = 0;
			s[1] = 0;
		}

		return s;
	}
}
