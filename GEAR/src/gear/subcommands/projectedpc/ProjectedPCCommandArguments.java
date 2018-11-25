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
			FileUtil.exists(tokens[0]+".bed");
			FileUtil.exists(tokens[0]+".bim");
			FileUtil.exists(tokens[0]+".fam");
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

	public void setInbred() {
		inbred = true;
	}

	public boolean isInbred() {
		return inbred;
	}

	public ArrayList<String> getBedFile() {
		return bedFiles;
	}

	public ArrayList<String> getAllBedFiles() {
		ArrayList<String> AllbedFiles = NewIt.newArrayList();
		AllbedFiles.add(getBFile());
		AllbedFiles.addAll(bedFiles);
		return AllbedFiles;
	}

	private String batch;
	private ArrayList<String> bedFiles = NewIt.newArrayList();

	private int ev = 2;
	private boolean inbred = false;
}
