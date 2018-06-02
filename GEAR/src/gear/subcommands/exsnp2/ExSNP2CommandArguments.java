package gear.subcommands.exsnp2;

import java.util.ArrayList;
import gear.subcommands.CommandArguments;
import gear.util.BufferedReader;
import gear.util.FileUtil;
import gear.util.NewIt;

public class ExSNP2CommandArguments extends CommandArguments {

	public void setBatch(String batch) {
		FileUtil.exists(batch);
		BatchFile = batch;

		BufferedReader reader = BufferedReader.openTextFile(batch, "Extract SNP Batch");

		String[] tokens = null;
		while ((tokens = reader.readTokens()) != null) {
			String sbed = tokens[0] + ".bed";
			FileUtil.exists(sbed);
			String sbim = tokens[0] + ".bim";
			FileUtil.exists(sbim);
			String sfam = tokens[0] + ".fam";
			FileUtil.exists(sfam);
			bFile.add(tokens[0]);
		}
	}

	public String getBatchFile() {
		return BatchFile;
	}

	public void setBFiles(String[] bf) {
		for (int i = 0; i < bf.length; i++) {
			String sbed = bf[i] + ".bed";
			FileUtil.exists(sbed);
			String sbim = bf[i] + ".bim";
			FileUtil.exists(sbim);
			String sfam = bf[i] + ".fam";
			FileUtil.exists(sfam);
			bFile.add(bf[i]);
		}
	}

	public ArrayList<String> getBFiles() {
		return bFile;
	}

	public String BatchFile = null;
	public ArrayList<String> bFile = NewIt.newArrayList();

}
