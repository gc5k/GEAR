package gear.subcommands.projectedpc;

import java.util.ArrayList;

import gear.subcommands.CommandArguments;
import gear.util.BufferedReader;
import gear.util.FileUtil;
import gear.util.NewIt;

public class ProjectedPCCommandArguments extends CommandArguments {

	public void setBatch(String batch) {
		FileUtil.exists(batch);
		this.batch = batch;

		BufferedReader reader = BufferedReader.openTextFile(this.batch, "Keep-cohort");

		String[] tokens = null;
		while ((tokens = reader.readTokens()) != null) {
			bedFiles.add(tokens[0]);
		}

	}

	public String getBatch() {
		return batch;
	}

	public void setEV(int ev) {
		if (ev < 1) {
			System.exit(1);
		}
		this.ev = ev;
	}

	public int getEV() {
		return ev;
	}

	private String batch;
	private ArrayList<String> bedFiles = NewIt.newArrayList();
	private int ev = 2;
}
