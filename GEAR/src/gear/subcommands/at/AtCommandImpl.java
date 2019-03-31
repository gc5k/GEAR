package gear.subcommands.at;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.PrintStream;

import gear.ConstValues;
import gear.family.GenoMatrix.GenotypeMatrix;
import gear.family.pedigree.file.SNP;
import gear.family.plink.PLINKParser;
import gear.qc.sampleqc.SampleFilter;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.FileUtil;
import gear.util.Logger;

public class AtCommandImpl extends CommandImpl {
	private AtCommandArguments atArgs = null;

	private SampleFilter sf;
	private GenotypeMatrix pGM;

	@Override
	public void execute(CommandArguments cmdArgs) {
		atArgs = (AtCommandArguments) cmdArgs;

		PLINKParser pp = PLINKParser.parse(this.atArgs);

		sf = new SampleFilter(pp.getPedigreeData());
		sf.qualification();
		pGM = new GenotypeMatrix(sf.getSample(), pp.getMapData(), cmdArgs);
		CalLDMat();
	}

	private void CalLDMat() {
		String ldOut = null;
		String sb_id = null;
		if (atArgs.isD()) {
			ldOut = new String(atArgs.getOutRoot() + ".d.gz");
			sb_id = new String(atArgs.getOutRoot() + ".d.snp");
			Logger.printUserLog("Calculating pairwise ld using d (covariance).");
		} else if (atArgs.isDPrime()) {
			ldOut = new String(atArgs.getOutRoot() + ".dprime.gz");
			sb_id = new String(atArgs.getOutRoot() + ".dprime.snp");
			Logger.printUserLog("Calculating pairwise ld using dprime (Lewontin's).");
		} else if (atArgs.isRsq()) {
			ldOut = new String(atArgs.getOutRoot() + ".rsq.gz");
			sb_id = new String(atArgs.getOutRoot() + ".rsq.snp");
			Logger.printUserLog("Calculating pairwise ld using rsq (squred correlation).");
		} else {
			ldOut = new String(atArgs.getOutRoot() + ".r.gz");
			sb_id = new String(atArgs.getOutRoot() + ".r.snp");
			Logger.printUserLog("Calculating pairwise ld using r (correlation).");
		}
		BufferedWriter mkGZ = FileUtil.ZipFileWriter(ldOut);

		long MkCnt = 0;
		for (int i = 0; i < pGM.getNumMarker(); i++) {
			SNP snp1 = pGM.getSNPList().get(i);
			for (int j = 0; j <= i; j++) {
				SNP snp2 = pGM.getSNPList().get(j);

				if (atArgs.isWindow()) {
					if (snp1.getChromosome().compareTo(snp2.getChromosome()) != 0)
						continue;
					if (Math.abs(snp1.getPosition() - snp2.getPosition()) > (atArgs.getWindow() * 1000000))
						continue;
				}

				MkCnt++;
				int cnt = 0;
				double sg1 = 0;
				double sg2 = 0;
				double ssg1 = 0;
				double ssg2 = 0;
				double crs = 0;
				float ldscore = 0;
				for (int k = 0; k < pGM.getNumIndivdial(); k++) {
					int g1 = pGM.getAdditiveScoreOnFirstAllele(k, i);
					int g2 = pGM.getAdditiveScoreOnFirstAllele(k, j);
					if (g1 == ConstValues.MISSING_GENOTYPE || g2 == ConstValues.MISSING_GENOTYPE) {
						continue;
					}
					sg1 += g1;
					sg2 += g2;
					ssg1 += g1 * g1;
					ssg2 += g2 * g2;
					crs += g1 * g2;
					cnt++;
				}
				if (cnt > 10) {
					double m1 = sg1 / (2 * cnt);
					double m2 = sg2 / (2 * cnt);
					double v1 = ssg1 / (4 * cnt) - m1 * m1;
					double v2 = ssg2 / (4 * cnt) - m2 * m2;
					double cv = crs / (4 * cnt) - m1 * m2;
					if (atArgs.isD()) {
						ldscore = (float) cv;
					} else if (atArgs.isDPrime()) {
						if (cv > 0) {
							ldscore = (float) (cv / Math.min(m1 * (1 - m2), m2 * (1 - m1)));
						} else {
							ldscore = (float) (cv / Math.min(m1 * m2, (1 - m2) * (1 - m1)));
						}
					} else if (atArgs.isRsq()) {
						ldscore = (float) (cv / Math.sqrt(v1 * v2));
						ldscore *= ldscore;
					} else {
						ldscore = (float) (cv / Math.sqrt(v1 * v2));
					}
				}
				try {
					mkGZ.append((i + 1) + "\t" + (j + 1) + "\t" + cnt + "\t" + ldscore + "\n");
				} catch (IOException e) {
					Logger.handleException(e,
							"error in writing '" + ldOut + "' for markers" + (i + 1) + " " + (j + 1) + ".");
				}
			}
		}
		try {
			mkGZ.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		Logger.printUserLog("Saved " + MkCnt + " ld scores for pairs of markers in '" + ldOut + "'.");

		PrintStream ld_mk = FileUtil.CreatePrintStream(sb_id);

		for (int i = 0; i < pGM.getNumMarker(); i++) {
			SNP snp = pGM.getSNPList().get(i);
			ld_mk.println(snp.getName() + " " + snp.getChromosome() + " " + snp.getDistance() + " " + snp.getPosition()
					+ " " + snp.getFirstAllele() + " " + snp.getSecAllele());
		}
		ld_mk.close();
		Logger.printUserLog("Saved " + pGM.getNumMarker() + " SNP information in '" + sb_id + "'.");
	}

}
