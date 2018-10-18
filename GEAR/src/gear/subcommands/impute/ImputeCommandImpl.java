package gear.subcommands.impute;

import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import gear.ConstValues;
import gear.data.InputDataSet2;
import gear.family.GenoMatrix.GenotypeMatrix;
import gear.family.pedigree.PersonIndex;
import gear.family.plink.PLINKBinaryParser;
import gear.family.plink.PLINKParser;
import gear.family.qc.rowqc.SampleFilter;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.subcommands.ibd.IBDCommandArguments;
import gear.util.Logger;
import gear.util.pop.PopStat;

public class ImputeCommandImpl extends CommandImpl {

	private ImputeCommandArguments imputeArgs;
	private InputDataSet2 data = new InputDataSet2();
	private SampleFilter sf;
	private GenotypeMatrix gm;

	@Override
	public void execute(CommandArguments cmdArgs) {
		imputeArgs = (ImputeCommandArguments) cmdArgs;
		PLINKParser pp = PLINKBinaryParser.parse(cmdArgs);

		data.addFile(imputeArgs.getFam()); // geno
		if (imputeArgs.isKeepFile())
			data.addFile(imputeArgs.getKeepFile()); // keep
		data.LineUpFiles();

		sf = new SampleFilter(pp.getPedigreeData(), data.getMatchSubjetList());
		gm = new GenotypeMatrix(sf.getSample(), pp.getMapData(), cmdArgs);

		naiveImpute();
		writeBfile();
	}

	private void writeBfile() {
		DataOutputStream bedout = null;
		PrintWriter fam = null;
		PrintWriter bim = null;
		try {
			bedout = new DataOutputStream(new FileOutputStream(imputeArgs.getOutRoot() + ".bed"));
			fam = new PrintWriter(new BufferedWriter(new FileWriter(imputeArgs.getOutRoot() + ".fam")));
			bim = new PrintWriter(new BufferedWriter(new FileWriter(imputeArgs.getOutRoot() + ".bim")));
		} catch (IOException e) {
			Logger.handleException(e, "An I/O exception occurred when creating the .bed, .fam and .bim files.");
		}

		for (int i = 0; i < gm.getNumIndivdial(); i++) {
			PersonIndex pi = sf.getSample().get(i);
			fam.print(pi.getFamilyID() + "\t");
			fam.print(pi.getIndividualID() + "\t");
			fam.print(pi.getPerson().getDadID() + "\t");
			
			fam.print(pi.getPerson().getMomID() + "\t");
			fam.print(pi.getPerson().getGender() + "\t");
			fam.println(pi.getPerson().getAffectedStatus());
		}
		fam.close();

		try {
			bedout.writeByte(ConstValues.PLINK_BED_BYTE1);
			bedout.writeByte(ConstValues.PLINK_BED_BYTE2);
			bedout.writeByte(ConstValues.PLINK_BED_BYTE3);
			for (int i = 0; i < gm.getNumMarker(); i++) {
				byte gbyte = 0;
				int idx = 0;
				for (int j = 0; j < gm.getNumIndivdial(); j++) {
					int g = (int) gm.getAdditiveScore(j, i);
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

					if (j != (gm.getNumIndivdial() - 1)) {
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

		for (int i = 0; i < gm.getNumMarker(); i++) {
			bim.print(gm.getSNPList().get(i).getChromosome() + " ");
			bim.print(gm.getSNPList().get(i).getName() + " ");
			bim.print(gm.getSNPList().get(i).getDistance() + " ");
			bim.print(gm.getSNPList().get(i).getPosition() + " ");
			bim.println(gm.getSNPList().get(i).getFirstAllele() + " " + gm.getSNPList().get(i).getSecAllele());
		}

		bim.close();

	}

	private void naiveImpute() {
		PopStat.Imputation(gm, imputeArgs.isInbred(), imputeArgs.getSeed());
	}

}
