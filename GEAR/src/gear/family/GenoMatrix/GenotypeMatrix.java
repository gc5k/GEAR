package gear.family.GenoMatrix;

import java.util.ArrayList;

import gear.CmdArgs;
import gear.ConstValues;
import gear.data.Person;
import gear.family.pedigree.PersonIndex;
import gear.family.pedigree.file.MapFile;
import gear.family.pedigree.file.SNP;
import gear.subcommands.CommandArguments;
import gear.util.Logger;
import gear.util.NewIt;
import gear.util.pop.PopStat;

public class GenotypeMatrix {
	private int[][] genotypeMat;
	private final int shift = 4;
	private int numMarker = 0;
	private ArrayList<PersonIndex> pidx;
	private ArrayList<SNP> snpList;
	private ArrayList<Integer> QCedSnpIndex = NewIt.newArrayList();
	private float[][] allelefreq;
	private float[] allelevar;

	public GenotypeMatrix(ArrayList<PersonIndex> pi, MapFile mapF) {
		QCedSnpIndex = NewIt.newArrayList();
		for (int i = 0; i < mapF.getMarkerNumber(); i++) {
			QCedSnpIndex.add(i);
		}
		snpList = mapF.getMarkerList();
		pidx = pi;
		genotypeMat = new int[pidx.size()][];

		initial();
		
		if (this.getNumMarker() == 0) {
			Logger.printUserLog("No SNPs were remained for analysis. GEAR quit.");
			System.exit(1);
		}
	}

	public GenotypeMatrix(ArrayList<PersonIndex> pi, MapFile mapF, CommandArguments cmdArgs) {
		QCedSnpIndex = NewIt.newArrayList();
		for (int i = 0; i < mapF.getMarkerNumber(); i++) {
			QCedSnpIndex.add(i);
		}
		snpList = mapF.getMarkerList();
		pidx = pi;
		genotypeMat = new int[pidx.size()][];

		initial(cmdArgs);

		snpList = NewIt.newArrayList();
		for (int i = 0; i < QCedSnpIndex.size(); i++) {
			snpList.add(mapF.getMarkerList().get(QCedSnpIndex.get(i)));
		}

		if (this.getNumMarker() == 0) {
			Logger.printUserLog("No SNPs were remained for analysis. GEAR quit.");
			System.exit(1);
		}

	}

	protected void initial(CommandArguments cmdArgs) {
		int c1 = 0;
		for (PersonIndex pi : pidx) {
			if (!pi.isPseudo()) {
				genotypeMat[c1++] = pi.getPerson().getAlleleArray();
			}
		}
		numMarker = pidx.get(0).getPerson().getNumMarkers();
		
		allelefreq = PopStat.calAlleleFrequencyFloat(this);
		allelevar = PopStat.calGenoVarianceFloat(this);

		if (cmdArgs.isMAF() || cmdArgs.isMaxMAF() || cmdArgs.isGENO() || cmdArgs.isZeroVar() || cmdArgs.isMAFRange()) {
			Logger.printUserLog("");
			Logger.printUserLog("Quality control for SNPs...");

			QCedSnpIndex.clear();
			int mafFail = 0;
			int maxmafFail = 0;
			int genoFail = 0;
			int zeroVarFail = 0;
			int mafRangeFail = 0;
			for (int i = 0; i < numMarker; i++) {
				double m = allelefreq[i][0] < 0.5 ? allelefreq[i][0] : (1 - allelefreq[i][0]);
				if (cmdArgs.isMAFRange()) {
					double[][] mafR = cmdArgs.getMAFRange();
					boolean mafFlag = false;
					for (int j = 0; j < mafR.length; j++) {
						if (m >= mafR[j][0] && m <= mafR[j][1]) {
							mafFlag = true;
						}
					}
					if (!mafFlag) {
						mafRangeFail++;
						continue;
					}
				}

				if (cmdArgs.isMAF() && m < cmdArgs.getMAF()) {
					mafFail++;
					continue;
				}
				if (cmdArgs.isMaxMAF() && m > cmdArgs.getMaxMAF()) {
					maxmafFail++;
					continue;
				}
				if (cmdArgs.isGENO() && allelefreq[i][2] > cmdArgs.getGENO()) {
					genoFail++;
					continue;
				}
				if (cmdArgs.isZeroVar() && allelevar[i] == 0) {
					zeroVarFail++;
					continue;
				}
				QCedSnpIndex.add(i);
			}

			if (cmdArgs.isMAFRange() && mafRangeFail > 0) {
				Logger.printUserLog(mafRangeFail + " SNPs were excluded because not within maf-range.");				
			}
			if (cmdArgs.isMAF() && mafFail > 0) {
				Logger.printUserLog(mafFail + " SNPs were excluded because did not pass maf threshold " + cmdArgs.getMAF() + ".");
			}
			if (cmdArgs.isMaxMAF() && maxmafFail > 0) {
				Logger.printUserLog(maxmafFail + " SNPs were excluded because did not pass max-maf threshold " + cmdArgs.getMaxMAF() + ".");
			}
			if (cmdArgs.isGENO() && genoFail > 0) {
				Logger.printUserLog(genoFail + " SNPs were excluded because did not pass geno missing threshold " + cmdArgs.getGENO() + ".");
			}
			if (cmdArgs.isZeroVar() && zeroVarFail > 0) {
				Logger.printUserLog(zeroVarFail + " SNPs were excluded because did not pass zero variance threshold.");
			}

			Logger.printUserLog(QCedSnpIndex.size() + " SNPs were remained for analysis.");
		}
	}

	protected void initial() {
		int c1 = 0;
		for (PersonIndex pi : pidx) {
			if (!pi.isPseudo()) {
				genotypeMat[c1++] = pi.getPerson().getAlleleArray();
			}
		}
		numMarker = pidx.get(0).getPerson().getNumMarkers();
		
		Logger.printUserLog("Calculating allele frequencies and variance...");

		allelefreq = PopStat.calAlleleFrequencyFloat(this);
		allelevar = PopStat.calGenoVarianceFloat(this);
	}

	public int getNumIndivdial() {
		return getGRow();
	}

	public int getNumMarker() {
		return QCedSnpIndex.size();
	}

	public void setAdditiveScore(int idx, int i, int v) {
		int posByte = QCedSnpIndex.get(i) >> shift;
		int posBit = (QCedSnpIndex.get(i) & 0xf) << 1;
		genotypeMat[idx][posByte] &= ~(3 << posBit);
		genotypeMat[idx][posByte] |= (v & 3) << posBit;
	}

	public boolean isNA(int idx, int i) {
		// 0 homozygote 1/1
		// 1 heterozygosity 1/2
		// 2 homozygote 2/2
		// 3 missing
		int posByte = QCedSnpIndex.get(i) >> shift;
		int posBit = (QCedSnpIndex.get(i) & 0xf) << 1;
		int g = (genotypeMat[idx][posByte] >> (posBit)) & 3;
		return g == 3? true:false;
		// if (g == 1) {// 01
		// return 3;
		// } else if (g == 0) {
		// return 0;
		// } else if (g == 2) {
		// return 1;
		// } else {
		// return 2;
		// }

	}

	public int getAdditiveScore(int idx, int i) {
		// 0 homozygote 1/1
		// 1 heterozygosity 1/2
		// 2 homozygote 2/2
		// 3 missing
		int posByte = QCedSnpIndex.get(i) >> shift;
		int posBit = (QCedSnpIndex.get(i) & 0xf) << 1;
		int g = (genotypeMat[idx][posByte] >> (posBit)) & 3;
		return g;
		// if (g == 1) {// 01
		// return 3;
		// } else if (g == 0) {
		// return 0;
		// } else if (g == 2) {
		// return 1;
		// } else {
		// return 2;
		// }
	}

	public int getAdditiveScoreOnFirstAllele(int idx, int i) {
		// 0 homozygote 1/1
		// 1 heterozygosity 1/2
		// 2 homozygote 2/2
		// 3 missing
		int posByte = QCedSnpIndex.get(i) >> shift;
		int posBit = (QCedSnpIndex.get(i) & 0xf) << 1;
		int g = (genotypeMat[idx][posByte] >> (posBit)) & 3;
		return g == ConstValues.MISSING_GENOTYPE ? g : (2 - g);
		// if (g == 1) {// 01
		// return 3;
		// } else if (g == 0) {
		// return 0;
		// } else if (g == 2) {
		// return 1;
		// } else {
		// return 2;
		// }
	}

	public byte getOriginalGenotypeScore(int ind, int i) {// this only for write back to bed file purpose
		int posByte = QCedSnpIndex.get(i) >> shift;
		int posBit = (QCedSnpIndex.get(i) & 0xf) << 1;
		int g = (genotypeMat[ind][posByte] >> (posBit)) & 3;
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
		return (byte) g;
	}

	public int[] getBiAlleleGenotype(int idx, int i) {
		int posByte = QCedSnpIndex.get(i) >> shift;
		int posBit = (QCedSnpIndex.get(i) & 0xf) << 1;
		int g = (genotypeMat[idx][posByte] >> posBit) & 3;
		int[] b = { 0, 0 };
		switch (g) {
		case 0:
			b[0] = 0;
			b[1] = 0;
			break;
		case 1:
			b[0] = 0;
			b[1] = 1;
			break;
		case 2:
			b[0] = 1;
			b[1] = 1;
			break;
		default:
			b[0] = Person.MissingAlleleCode;
			b[1] = Person.MissingAlleleCode;
			break;
		}
		return b;
	}

	public String getGenotypeScoreString(int idx, int i) {
		int posByte = QCedSnpIndex.get(i) >> shift;
		int posBit = (QCedSnpIndex.get(i) & 0xf) << 1;
		int g = (genotypeMat[idx][posByte] >> (posBit)) & 3;
		if (g == 3) {// 01
			return CmdArgs.INSTANCE.missingGenotype;
		} else {
			return Integer.toString(g);
		}
	}

	public int getGRow() {
		return genotypeMat.length;
	}

	public int getGCol() {
		return QCedSnpIndex.size();
		// return genotypeMat[0].length;
	}

	public int[][] getG() {
		return genotypeMat;
	}

	public ArrayList<SNP> getSNPList() {
		return snpList;
	}

	public ArrayList<PersonIndex> getSample() {
		return pidx;
	}

	public double getAlleleVar(int i) {
		return allelevar[QCedSnpIndex.get(i)];
	}

	public float[] getQCedGenoVar() {
		float[] QCedAlleleVar = new float[QCedSnpIndex.size()];
		for(int i = 0; i < QCedSnpIndex.size(); i++) {
			QCedAlleleVar[i] = allelevar[QCedSnpIndex.get(i)];
		}
		return QCedAlleleVar;
	}

	public double getMAF(int i) {
		return allelefreq[QCedSnpIndex.get(i)][0] < 0.5 ? allelefreq[QCedSnpIndex.get(i)][0]:allelefreq[QCedSnpIndex.get(i)][1];
	}

	public double getA1Freq(int i) {
		return allelefreq[QCedSnpIndex.get(i)][0];
	}

	public double getA2Freq(int i) {
		return allelefreq[QCedSnpIndex.get(i)][1];
	}

	public float[][] getQCedAlleleFreq() {
		float[][] QCedAlleleFreq = new float[QCedSnpIndex.size()][3];
		for(int i = 0; i < QCedSnpIndex.size(); i++) {
			System.arraycopy(allelefreq[QCedSnpIndex.get(i)], 0, QCedAlleleFreq[i], 0, 3);
		}
		return QCedAlleleFreq;
	}

	public double getGenoMissing(int i) {
		return allelefreq[QCedSnpIndex.get(i)][2];
	}

	public float[] getQCedGenoMissing() {
		float[] QCedGenoMissing = new float[QCedSnpIndex.size()];
		for(int i = 0; i < QCedSnpIndex.size(); i++) {
			QCedGenoMissing[i] = allelefreq[i][2];
		}
		return QCedGenoMissing;
	}
}
