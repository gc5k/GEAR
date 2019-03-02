package gear.subcommands.exsnp2;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import gear.family.GenoMatrix.GenotypeMatrix;
import gear.family.pedigree.file.SNP;
import gear.family.plink.PLINKParser;
import gear.qc.sampleqc.SampleFilter;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.BufferedReader;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;

public class ExSNP2CommandImpl extends CommandImpl {

	private ExSNP2CommandArguments esArgs;

	private GenotypeMatrix pGM;
	private SampleFilter sf;
	private HashMap<String, Integer> snpCnt = NewIt.newHashMap();
	private ArrayList<SNP> snpInfo = NewIt.newArrayList();

	@Override
	public void execute(CommandArguments cmdArgs) {
		
		esArgs = (ExSNP2CommandArguments) cmdArgs;

		PLINKParser pp = PLINKParser.parse(esArgs);
		sf = new SampleFilter(pp.getPedigreeData(), cmdArgs);
		pGM = new GenotypeMatrix(sf.getSample(), pp.getMapData(), cmdArgs);

		for (int i = 0; i < pGM.getNumMarker(); i++) {
			SNP snp = pGM.getSNPList().get(i);
				snpInfo.add(snp);
			if (snpCnt.containsKey(snp.getName())) {
				Integer cnt = snpCnt.get(snp.getName());
				cnt++;
				snpCnt.put(snp.getName(), cnt);
			} else {
				snpCnt.put(snp.getName(), 1);
			}
		}
		Logger.printUserLog("Read " + pGM.getNumMarker() + " snps in '" +  esArgs.getBim()+"'.");

		ArrayList<String> bFile = esArgs.getBFiles();
		if (bFile.size() == 0) {
			Logger.printUserLog("no bfiles found. GEAR quit.");
			System.exit(0);
		}
		for (int i = 0; i < bFile.size(); i++) {
			int SNPcnt = 0;
			String bim = bFile.get(i) + ".bim";
			BufferedReader reader = BufferedReader.openTextFile(bim, "Extract SNP Batch");

			String[] tokens = null;
			while ((tokens = reader.readTokens()) != null) {
				if (tokens.length == 6) {
					if (snpCnt.containsKey(tokens[1])) {
						Integer cnt = snpCnt.get(tokens[1]);
						cnt++;
						snpCnt.put(tokens[1], cnt);
					} else {
						snpCnt.put(tokens[1], 1);
					}
				}
				SNPcnt++;
			}
			Logger.printUserLog("Read " + SNPcnt + " snps in '" + bFile.get(i) + ".bim'.");
		}

		ArrayList<SNP> finalSNP = NewIt.newArrayList();

		for (int i = 0; i < snpInfo.size(); i++) {
			SNP snp = snpInfo.get(i);
			String key = snp.getName();
			Integer cnt = snpCnt.get(key);
			if (cnt.intValue() == (bFile.size() + 1)) {
				finalSNP.add(snp);
			}
		}

		Collections.sort(finalSNP);
		Logger.printUserLog("Found " + finalSNP.size() + " consensus snps.");

		StringBuffer sb = new StringBuffer();
		sb.append(esArgs.getOutRoot());
		sb.append(".snp");
		PrintStream ps = FileUtil.CreatePrintStream(sb.toString());

		for (int i = 0; i < finalSNP.size(); i++) {
			ps.println(finalSNP.get(i).getName());
		}
		ps.close();
		Logger.printUserLog("Save results in '" + esArgs.getOutRoot() + ".snp'.");
	}

}
