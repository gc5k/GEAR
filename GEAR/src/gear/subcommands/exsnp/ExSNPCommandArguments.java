package gear.subcommands.exsnp;

import java.util.ArrayList;
import gear.subcommands.CommandArguments;
import gear.util.BufferedReader;
import gear.util.FileUtil;
import gear.util.NewIt;

public class ExSNPCommandArguments extends CommandArguments {

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

	public void addBfile(String file) {
		String sbed = file + ".bed";
		FileUtil.exists(sbed);
		String sbim = file + ".bim";
		FileUtil.exists(sbim);
		String sfam = file + ".fam";
		FileUtil.exists(sfam);
		bFile.add(file);
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
