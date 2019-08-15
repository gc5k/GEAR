package gear.subcommands.simulationfam;

import gear.ConstValues;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.Sample;
import gear.util.pop.PopStat;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistributionImpl;
import org.apache.commons.math.random.RandomDataImpl;
import org.apache.commons.math.stat.StatUtils;

public final class SimuFamilyCommandImpl extends CommandImpl {
	@Override
	public void execute(CommandArguments cmdArgs) {
		famArgs = (SimuFamilyCommandArguments) cmdArgs;

		init();
		try {
			ibdF = new PrintWriter(new BufferedWriter(new FileWriter(famArgs.getOutRoot() + ".ibdo")));
			ibdSibF = new PrintWriter(new BufferedWriter(new FileWriter(famArgs.getOutRoot() + ".hap")));
		} catch (IOException e) {
			Logger.handleException(e, "An I/O exception occurred when creating the .ibdo file.");
			Logger.handleException(e, "An I/O exception occurred when creating the .ibd file.");
		}

		for (int i = 0; i < famArgs.getNumberOfFamilies(); i++) {
			generateNuclearFamily(NKid[i], NAffKid[i], i);
		}
		ibdF.close();
		ibdSibF.close();

		writeIBD();
		writePheno();
		writeEffFile();

		if (famArgs.getMakeBed()) {
			writeBFile();
		} else {
			writeFile();
		}
	}

	private void init() {
		rnd = new RandomDataImpl();
		rnd.reSeed(famArgs.getSeed());

		M = famArgs.getNumberOfMarkers();
		nullM = famArgs.getNullMarker();
		NKid = new int[famArgs.getNumberOfFamilies()];
		Arrays.fill(NKid, 2);
		NAffKid = new int[famArgs.getNumberOfFamilies()];
		Arrays.fill(NAffKid, 1);

		hsq = famArgs.getHsq();
		getEffect();
		getDomEffect();
		getMAF();
		getDPrime();
		calLD();
		getRec();

		gm = new int[famArgs.getNumberOfFamilies() * famSize][M];
		phe = new double[famArgs.getNumberOfFamilies() * famSize];
		BV = new double[famArgs.getNumberOfFamilies() * famSize];
		IBD = new double[freq.length][famArgs.getNumberOfFamilies()];
		IBDp = new double[freq.length][famArgs.getNumberOfFamilies()];
		IBDm = new double[freq.length][famArgs.getNumberOfFamilies()];
		IBD2 = new double[freq.length][famArgs.getNumberOfFamilies()];
	}

	private void getRec() {
		recSex = new double[M][2];
		if (famArgs.isPlainRec()) {
			for (int i = 0; i < recSex.length; i++) {
				recSex[i][0] = famArgs.getRec();
				recSex[i][1] = famArgs.getRec();
			}
		} else if (famArgs.isSexRec()) {
			double[] rs = famArgs.getRecSex();
			for (int i = 0; i < recSex.length; i++) {
				recSex[i][0] = rs[0];
				recSex[i][1] = rs[1];
			}
		} else if (famArgs.isRandRec()) {
			for (int i = 0; i < recSex.length; i++) {
				recSex[i][0] = recSex[i][1] = rnd.nextUniform(0.01, 0.5);
			}
		}
		recSex[0][0] = 0.5;
		recSex[0][1] = 0.5;
	}

	private void calLD() {
		LD = PopStat.CalcLDfromDPrime(freq, DPrime);
	}

	private void getDPrime() {
		DPrime = new double[M - 1];
		if (famArgs.isPlainLD()) {
			Arrays.fill(DPrime, famArgs.getLD());
		} else if (famArgs.isRandLD()) {
			for (int i = 0; i < DPrime.length; i++) {
				DPrime[i] = rnd.nextUniform(-1, 1);
			}
		}
	}

	private void getMAF() {
		freq = new double[M];
		if (famArgs.isPlainMAF()) {
			Arrays.fill(freq, famArgs.getMAF());
		} else if (famArgs.isUnifMAF()) {
			for (int i = 0; i < freq.length; i++) {
				freq[i] = rnd.nextUniform(0.01, 0.5);
			}
		} else if (famArgs.isFreqFile()) {
			BufferedReader reader = FileUtil.FileOpen(famArgs.getFreqFile());
			int c = 0;
			String line = null;
			try {
				while ((line = reader.readLine()) != null) {
					if (c >= M) {
						Logger.printUserLog(
								"Have already read " + M + " allelic frequencies.  Ignore the rest of the content in '"
										+ famArgs.getFreqFile() + "'.");
						break;
					}
					line.trim();
					String[] l = line.split(ConstValues.WHITESPACE_DELIMITER);
					if (l.length < 1)
						continue;
					freq[c++] = Double.parseDouble(l[0]);
				}
				reader.close();
			} catch (IOException e) {
				Logger.handleException(e,
						"An exception occurred when reading the frequency file '" + famArgs.getFreqFile() + "'.");
			}
		}
	}

	private void getEffect() {
		int[] idx = Sample.SampleIndex(0, M - 1, M - nullM);
		Arrays.sort(idx);

		effect = new double[M];
		if (famArgs.isPolyEffectFile()) {
			BufferedReader reader = FileUtil.FileOpen(famArgs.getPolyEffectFile());

			int c = 0;
			String line = null;
			try {
				while ((line = reader.readLine()) != null) {
					if (c >= M) {
						Logger.printUserLog(
								"Have already read " + M + " allelic effects. Ignore the rest of the content in '"
										+ famArgs.getPolyEffectFile() + "'.");
						break;
					}

					line.trim();
					String[] l = line.split(ConstValues.WHITESPACE_DELIMITER);
					if (l.length < 1)
						continue;
					effect[c++] = Double.parseDouble(l[0]);
				}
				reader.close();
			} catch (IOException e) {
				Logger.handleException(e,
						"An exception occurred when reading the effect file '" + famArgs.getPolyEffectFile() + "'.");
			}
		} else {
			for (int i = 0; i < idx.length; i++) {
				effect[idx[i]] = rnd.nextGaussian(0, 1);
			}
		}

		if (hsq == 0) {
			for (int i = 0; i < effect.length; i++) {
				effect[i] = 0;
			}
		}
	}

	private void getDomEffect() {
		deffect = new double[M];
		if (famArgs.getHsqDom() == 0) {
			Arrays.fill(deffect, 0);
			return;
		}

		Sample.setSeed(famArgs.getSeed() + 1);

		int[] idx = Sample.SampleIndex(0, M - 1, M - nullM);
		Arrays.sort(idx);

		if (famArgs.isPlainDomEffect()) {
			for (int i = 0; i < idx.length; i++)
				deffect[idx[i]] = famArgs.getPolyDomEffect();
		} else if (famArgs.isPolyDomEffect()) {
			double t = 0;
			double adj = 1;
			if (famArgs.isPolyDomEffect()) {
				if (famArgs.getHsq() > 0) {
					t = famArgs.getHsq() / famArgs.getHsqDom();
					if (t < 2) {
						Logger.printUserError("Impossible heritability: hsq_add = " + famArgs.getHsq() + " hsq_dom = "
								+ famArgs.getHsqDom());
						System.exit(0);
					} else {
						adj = Math.sqrt(5 / (2 * t - 1));
					}
				}
			}

			for (int i = 0; i < idx.length; i++) {
				deffect[idx[i]] = rnd.nextGaussian(0, adj);
			}
		} else if (famArgs.isPolyDomEffectSort()) {
			NormalDistributionImpl ndImpl = new NormalDistributionImpl();
			ndImpl.reseedRandomGenerator(famArgs.getSeed() + 1);
			for (int i = 0; i < idx.length; i++) {
				try {
					deffect[idx[i]] = ndImpl.inverseCumulativeProbability((i + 0.5) / (idx.length));
				} catch (MathException e) {
					e.printStackTrace();
				}
			}
		} else if (famArgs.isPolyDomEffectFile()) {
			BufferedReader reader = FileUtil.FileOpen(famArgs.getPolyDomEffectFile());
			int c = 0;
			String line = null;
			try {
				while ((line = reader.readLine()) != null) {
					if (c >= M) {
						Logger.printUserLog(
								"Have already read " + M + " allelic effects.  Ignore the rest of the content in '"
										+ famArgs.getFreqFile() + "'.");
						break;
					}

					line.trim();
					String[] l = line.split(ConstValues.WHITESPACE_DELIMITER);
					if (l.length < 1)
						continue;
					deffect[c++] = Double.parseDouble(l[0]);
				}
				reader.close();
			} catch (IOException e) {
				Logger.handleException(e, "An exception occurred when reading the frequency file '"
						+ famArgs.getPolyDomEffectFile() + "'.");
			}
		}
	}

	private void generateNuclearFamily(int nkid, int affKid, int famIdx) {
		int[][] p = sampleChromosome(famIdx, 0);
		int[][] m = sampleChromosome(famIdx, 1);

		int[][] rc1 = generateBaby(p, m, famIdx, 2);
		for (int j = 0; j < 2; j++) {
			for (int k = 0; k < rc1.length; k++) {
				ibdF.print(rc1[k][j] + " ");
			}
			ibdF.println();
		}

		int[][] rc2 = generateBaby(p, m, famIdx, 3);
		for (int j = 0; j < 2; j++) {
			for (int k = 0; k < rc2.length; k++) {
				ibdF.print(rc2[k][j] + " ");
			}
			ibdF.println();
		}

		double[] ibd = new double[freq.length];
		for (int i = 0; i < freq.length; i++) {
			ibd[i] = 2 - Math.abs(rc1[i][0] - rc2[i][0]) - Math.abs(rc1[i][1] - rc2[i][1]);
			IBD[i][famIdx] = ibd[i] / 2;
			IBDp[i][famIdx] = 1 - Math.abs(rc1[i][0] - rc2[i][0]);
			IBDm[i][famIdx] = 1 - Math.abs(rc1[i][1] - rc2[i][1]);
			IBD2[i][famIdx] = (1 - Math.abs(rc1[i][0] - rc2[i][0])) * (1 - Math.abs(rc1[i][1] - rc2[i][1]));
			ibdSibF.print(ibd[i] + " ");
		}
		ibdSibF.println();
	}

	private int[][] sampleChromosome(int famIdx, int shift) {
		int[][] v = new int[freq.length][2];
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < freq.length; j++) {
				double r = rnd.nextUniform(0, 1);
				if (j == 0) {
					v[j][i] = r < freq[j] ? 0 : 1;
				} else {
					double d = rnd.nextUniform(0, 1);
					int a = (int) v[j - 1][i];
					double f1 = a == 0 ? freq[j - 1] : (1 - freq[j - 1]);
					double f2 = a == 0 ? freq[j] : (1 - freq[j]);
					v[j][i] = d < (f1 * f2 + LD[j - 1]) / f1 ? v[j - 1][i] : (1 - v[j - 1][i]);
				}
			}
		}

		for (int i = 0; i < freq.length; i++) {
			gm[famIdx * famSize + shift][i] = v[i][0] + v[i][1];
		}

		for (int i = 0; i < freq.length; i++) {
			phe[famIdx * famSize + shift] += gm[famIdx * famSize + shift][i] * effect[i];
			phe[famIdx * famSize + shift] += gm[famIdx * famSize + shift][i] == 1 ? deffect[i] : 0;
		}
		BV[famIdx * famSize + shift] = phe[famIdx * famSize + shift];

		return v;
	}

	private int[][] generateBaby(int[][] p, int[][] m, int famIdx, int shift) {
		int[][] v = new int[freq.length][2];
		int[][] rc = new int[freq.length][2];

		for (int i = 0; i < 2; i++) {
			int[][] chr = i == 0 ? p : m;
			int idx = 1;
			try {
				idx = rnd.nextBinomial(1, recSex[0][i]);
			} catch (MathException e) {
				e.printStackTrace();
			}

			for (int j = 0; j < freq.length; j++) {
				double r = rnd.nextUniform(0, 1);
				idx = r < recSex[j][i] ? 1 - idx : idx;
				rc[j][i] = idx;
				v[j][i] = chr[j][idx];
			}
		}

		for (int i = 0; i < freq.length; i++) {
			gm[famIdx * famSize + shift][i] = v[i][0] + v[i][1];
		}

		for (int i = 0; i < freq.length; i++) {
			phe[famIdx * famSize + shift] += gm[famIdx * famSize + shift][i] * effect[i];
			phe[famIdx * famSize + shift] += gm[famIdx * famSize + shift][i] == 1 ? deffect[i] : 0;
		}
		BV[famIdx * famSize + shift] = phe[famIdx * famSize + shift];

		return rc;
	}

	public void writeIBD() {
		StringBuffer sb = new StringBuffer();
		sb.append(famArgs.getOutRoot() + ".ibd.gz");
		BufferedWriter ibdGZFile = FileUtil.ZipFileWriter(sb.toString());

		StringBuffer sbp = new StringBuffer();
		sbp.append(famArgs.getOutRoot() + ".p.ibd.gz");
		BufferedWriter ibdpGZFile = FileUtil.ZipFileWriter(sbp.toString());

		StringBuffer sbm = new StringBuffer();
		sbm.append(famArgs.getOutRoot() + ".m.ibd.gz");
		BufferedWriter ibdmGZFile = FileUtil.ZipFileWriter(sbm.toString());

		StringBuffer sbDom = new StringBuffer();
		sbDom.append(famArgs.getOutRoot() + ".dom.ibd.gz");
		BufferedWriter ibdDomGZFile = FileUtil.ZipFileWriter(sbDom.toString());

		// for (int i = 0; i < freq.length; i++)
		// {
		// IBD[i] = StatUtils.normalize(IBD[i]);
		// IBDp[i] = StatUtils.normalize(IBDp[0]);
		// IBDm[i] = StatUtils.normalize(IBDm[0]);
		// IBD2[i] = StatUtils.normalize(IBD2[i]);
		// }

		for (int i = 0; i < famArgs.getNumberOfFamilies(); i++) {
			double ibdA = 0;
			double ibdAp = 0;
			double ibdAm = 0;
			double ibdD = 0;
			for (int j = 0; j < freq.length; j++) {
				ibdA += IBD[j][i];
				ibdAp += IBDp[j][i];
				ibdAm += IBDm[j][i];
				ibdD += IBD2[j][i];
			}
			ibdA /= freq.length;
			ibdAm /= freq.length;
			ibdAp /= freq.length;
			ibdD /= freq.length;
			try {
				ibdGZFile.append((i + 1) * 4 + " " + ((i + 1) * 4 - 1) + " " + freq.length + " " + ibdA + "\n");
				ibdpGZFile.append((i + 1) * 4 + " " + ((i + 1) * 4 - 1) + " " + freq.length + " " + ibdAp + "\n");
				ibdmGZFile.append((i + 1) * 4 + " " + ((i + 1) * 4 - 1) + " " + freq.length + " " + ibdAm + "\n");
				ibdDomGZFile.append((i + 1) * 4 + " " + ((i + 1) * 4 - 1) + " " + freq.length + " " + ibdD + "\n");
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		try {
			ibdGZFile.close();
			ibdpGZFile.close();
			ibdmGZFile.close();
			ibdDomGZFile.close();
		} catch (IOException e1) {
			e1.printStackTrace();
		}

		PrintWriter famA = null;
		PrintWriter famAp = null;
		PrintWriter famAm = null;
		PrintWriter famD = null;

		try {
			famA = new PrintWriter(new BufferedWriter(new FileWriter(famArgs.getOutRoot() + ".ibd.id")));
			famAp = new PrintWriter(new BufferedWriter(new FileWriter(famArgs.getOutRoot() + ".p.ibd.id")));
			famAm = new PrintWriter(new BufferedWriter(new FileWriter(famArgs.getOutRoot() + ".m.ibd.id")));
			famD = new PrintWriter(new BufferedWriter(new FileWriter(famArgs.getOutRoot() + ".dom.ibd.id")));
		} catch (IOException e) {
			Logger.handleException(e, "An I/O exception occurred when creating the .ibd.id or .dom.ibd.id files.");
		}

		int fid = 0;
		for (int i = 1; i <= famArgs.getNumberOfFamilies(); i++) {
			fid = i * 10000;
			famA.println(fid + " " + (fid + 1));
			famA.println(fid + " " + (fid + 2));
			famA.println(fid + " " + (fid + 3));
			famA.println(fid + " " + (fid + 4));

			famAm.println(fid + " " + (fid + 1));
			famAm.println(fid + " " + (fid + 2));
			famAm.println(fid + " " + (fid + 3));
			famAm.println(fid + " " + (fid + 4));

			famAp.println(fid + " " + (fid + 1));
			famAp.println(fid + " " + (fid + 2));
			famAp.println(fid + " " + (fid + 3));
			famAp.println(fid + " " + (fid + 4));

			famD.println(fid + " " + (fid + 1));
			famD.println(fid + " " + (fid + 2));
			famD.println(fid + " " + (fid + 3));
			famD.println(fid + " " + (fid + 4));
		}
		famA.close();
		famAm.close();
		famAp.close();
		famD.close();
	}

	public void writePheno() {
		PrintWriter pheno = null;
		try {
			pheno = new PrintWriter(new BufferedWriter(new FileWriter(famArgs.getOutRoot() + ".phe")));
		} catch (IOException e) {
			Logger.handleException(e, "An I/O exception occurred when creating the .phe file.");
		}

		double[] p = new double[famArgs.getNumberOfFamilies() * 2];
		double[] po = new double[famArgs.getNumberOfFamilies()];
		int cn = 0;
		int cno = 0;

		for (int h = 0; h < famArgs.getNumberOfFamilies(); h++) {
			for (int j = 0; j < famSize; j++) {
				if (j < 2) {
					p[cn++] = phe[h * famSize + j];
				} else if (j == 2) {
					po[cno++] = phe[h * famSize + j];
				}
			}
		}

		double vg = StatUtils.variance(p);
		double vgo = StatUtils.variance(po);

		Logger.printUserLog("Genetic variation for founders is " + vg);
		Logger.printUserLog("Genetic variation for non-fouders is " + vgo);

		double H2 = famArgs.getHsq() + famArgs.getHsqDom();
		double ve = H2 == 0 ? 1 : vg * (1 - H2) / H2;

		for (int h = 0; h < famArgs.getNumberOfFamilies(); h++) {
			int fid = (h + 1) * 10000;

			for (int j = 0; j < famSize; j++) {
				int pid = fid + 1 + j;

				pheno.print(fid + " ");
				pheno.print(pid + " ");
				for (int r = 0; r < famArgs.getRep(); r++) {
					double pv = phe[h * famSize + j] + rnd.nextGaussian(0, Math.sqrt(ve));
					pheno.print(pv + " ");
				}
				pheno.println();
			}
		}

		pheno.close();
	}

	public void writeFile() {
		PrintWriter ped = null;
		PrintWriter map = null;
		PrintWriter breed = null;

		try {
			ped = new PrintWriter(new BufferedWriter(new FileWriter(famArgs.getOutRoot() + ".ped")));
			map = new PrintWriter(new BufferedWriter(new FileWriter(famArgs.getOutRoot() + ".map")));
			breed = new PrintWriter(new BufferedWriter(new FileWriter(famArgs.getOutRoot() + ".breed")));
		} catch (IOException e) {
			Logger.handleException(e, "An I/O exception occurred when creating the .ped and .map files.");
		}

		for (int h = 0; h < famArgs.getNumberOfFamilies(); h++) {
			int fid = (h + 1) * 10000;
			int pid = fid + 1;
			int mid = fid + 2;

			ped.print(fid + " ");
			ped.print(pid + " ");
			ped.print(0 + " ");
			ped.print(0 + " ");
			ped.print(1 + " ");
			ped.print(1 + " ");
			breed.println(fid + " " + pid + " " + BV[h * 4]);

			for (int i = 0; i < freq.length; i++) {
				StringBuilder sb = new StringBuilder();
				switch (gm[h * famSize][i]) {
				case 0:
					sb.append(A1 + " " + A1);
					break;
				case 1:
					sb.append(A1 + " " + A2);
					break;
				case 2:
					sb.append(A2 + " " + A2);
					break;
				default:
					break;
				}
				if (i == (freq.length - 1)) {
					ped.print(sb.toString());
				} else {
					ped.print(sb.toString() + " ");
				}
			}
			ped.print("\n");

			ped.print(fid + " ");
			ped.print(mid + " ");
			ped.print(0 + " ");
			ped.print(0 + " ");
			ped.print(2 + " ");
			ped.print(1 + " ");
			breed.println(fid + " " + mid + " " + BV[h * 4 + 1]);

			for (int i = 0; i < freq.length; i++) {
				StringBuilder sb = new StringBuilder();
				switch (gm[h * famSize + 1][i]) {
				case 0:
					sb.append(A1 + " " + A1);
					break;
				case 1:
					sb.append(A1 + " " + A2);
					break;
				case 2:
					sb.append(A2 + " " + A2);
					break;
				default:
					break;
				}
				if (i == (freq.length - 1)) {
					ped.print(sb.toString());
				} else {
					ped.print(sb.toString() + " ");
				}
			}
			ped.print("\n");

			for (int j = 0; j < 2; j++) {
				ped.print(fid + " ");
				ped.print((fid + 3 + j) + " ");
				ped.print(pid + " ");
				ped.print(mid + " ");
				breed.println(fid + " " + (fid + 3 + j) + " " + BV[h * 4 + 2 + j]);

				try {
					ped.print((rnd.nextBinomial(1, 0.5) + 1) + " ");
					ped.print((rnd.nextBinomial(1, 0.5) + 1) + " ");
				} catch (MathException e) {
					Logger.handleException(e, "Failed to generate the random values.");
				}

				for (int i = 0; i < freq.length; i++) {
					StringBuilder sb = new StringBuilder();
					switch (gm[h * famSize + 2 + j][i]) {
					case 0:
						sb.append(A1 + " " + A1);
						break;
					case 1:
						sb.append(A1 + " " + A2);
						break;
					case 2:
						sb.append(A2 + " " + A2);
						break;
					default:
						break;
					}
					if (i == (freq.length - 1)) {
						ped.print(sb.toString());
					} else {
						ped.print(sb.toString() + " ");
					}
				}
				ped.print("\n");
			}
		}

		for (int i = 0; i < freq.length; i++) {
			map.print(1 + " ");
			map.print("rs" + i + " ");
			map.print(i / (freq.length * 1.0) + " ");
			map.println(i * 100);
		}

		ped.close();
		map.close();
	}

	public void writeBFile() {
		DataOutputStream bedout = null;
		PrintWriter fam = null;
		PrintWriter bim = null;
		PrintWriter breed = null;
		try {
			bedout = new DataOutputStream(new FileOutputStream(famArgs.getOutRoot() + ".bed"));
			fam = new PrintWriter(new BufferedWriter(new FileWriter(famArgs.getOutRoot() + ".fam")));
			bim = new PrintWriter(new BufferedWriter(new FileWriter(famArgs.getOutRoot() + ".bim")));
			breed = new PrintWriter(new BufferedWriter(new FileWriter(famArgs.getOutRoot() + ".breed")));
		} catch (IOException e) {
			Logger.handleException(e, "An I/O exception occurred when creating the .bed, .fam and .bim files.");
		}

		int fid = 0;
		int fc = 0;
		int cn = 0;
		int pid = 0;
		int mid = 0;
		for (int i = 0; i < gm.length; i++) {
			if (i % 4 == 0) {
				fc++;
				fid = fc * 10000;
				pid = fid + 1;
				mid = fid + 2;
				cn = 0;
			}
			cn++;

			breed.println(fid + " " + (fid + cn) + " " + BV[i]);

			fam.print(fid + " ");
			if (cn <= 2) {
				fam.print((fid + cn) + " ");
				fam.print(0 + " ");
				fam.print(0 + " ");
			} else {
				fam.print((fid + cn) + " ");
				fam.print(pid + " ");
				fam.print(mid + " ");
			}
			try {
				if (cn <= 2) {
					fam.print(cn + " ");
				} else {
					fam.print((rnd.nextBinomial(1, 0.5) + 1) + " ");
				}
				fam.print((rnd.nextBinomial(1, 0.5) + 1) + "\n");

			} catch (MathException e) {
				Logger.handleException(e, "Failed to generate the random values.");
			}
		}

		try {
			bedout.writeByte(ConstValues.PLINK_BED_BYTE1);
			bedout.writeByte(ConstValues.PLINK_BED_BYTE2);
			bedout.writeByte(ConstValues.PLINK_BED_BYTE3);
			for (int i = 0; i < M; i++) {
				byte gbyte = 0;
				int idx = 0;
				for (int j = 0; j < gm.length; j++) {
					int g = (int) gm[j][i];
					switch (g) {
					case 0:
						g = 0;
						break;
					case 1:
						g = 2;
						break;
					case 2:
						g = 3;
						break;
					default:
						g = 1;
						break; // missing
					}

					g <<= 2 * idx;
					gbyte |= g;
					idx++;

					if (j != (gm.length - 1)) {
						if (idx == 4) {
							bedout.writeByte(gbyte);
							gbyte = 0;
							idx = 0;
						}
					} else {
						bedout.writeByte(gbyte);
					}
				}
			}
			bedout.close();
		} catch (IOException e) {
			Logger.handleException(e, "An I/O exception occurred when writing the .bed file.");
		}

		for (int i = 0; i < M; i++) {
			bim.print(1 + " ");
			bim.print("rs" + i + " ");
			bim.print(i / (M * 1.0) + " ");
			bim.print(i * 100 + " ");
			bim.println(A1 + " " + A2);
		}

		bim.close();
		fam.close();
		breed.close();
	}

	private void writeEffFile() {
		PrintWriter eff = null;
		PrintWriter deff = null;
		try {
			eff = new PrintWriter(new BufferedWriter(new FileWriter(famArgs.getOutRoot() + ".rnd")));
			deff = new PrintWriter(new BufferedWriter(new FileWriter(famArgs.getOutRoot() + ".drnd")));
		} catch (IOException e) {
			Logger.handleException(e, "An exception occurred when writing files.");
		}

		for (int i = 0; i < M; i++) {
			eff.println("rs" + i + " " + A1 + " " + effect[i]);
			deff.println("rs" + i + " " + deffect[i]);
		}
		eff.close();
		deff.close();
	}

	private SimuFamilyCommandArguments famArgs;

	private RandomDataImpl rnd;
	private int[] NKid = null;
	private int[] NAffKid = null;
	private double[] LD = null;
	private double[][] recSex = null;
	private double[] freq = null;
	private double[] DPrime = null;

	private int[][] gm = null;
	private double[] BV = null;
	private final int famSize = 4;
	private double[] phe = null;

	private int M;
	private int nullM;

	private double hsq;
	private double[] effect;
	private double[] deffect;

	PrintWriter ibdF = null;
	PrintWriter ibdSibF = null;

	private double[][] IBD = null;
	private double[][] IBDp = null;
	private double[][] IBDm = null;

	private double[][] IBD2 = null;

	private String A1 = "A";
	private String A2 = "C";

}
