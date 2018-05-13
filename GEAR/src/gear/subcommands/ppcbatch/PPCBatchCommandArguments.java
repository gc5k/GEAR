package gear.subcommands.ppcbatch;

import java.util.ArrayList;

import gear.subcommands.CommandArguments;
import gear.subcommands.exsnp.ExSNPCommandArguments;
import gear.subcommands.profile.ProfileCommandArguments;
import gear.util.BufferedReader;
import gear.util.FileUtil;
import gear.util.NewIt;

public class PPCBatchCommandArguments extends CommandArguments {

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

	public ArrayList<String> getBedFiles() {
		return bedFiles;
	}

	public void setHeader() {
		hasHeader = false;
	}

	public boolean getHeader() {
		return hasHeader;
	}

	public void setScoreFile(String sFile) {
		FileUtil.exists(sFile);
		scoreFile = sFile;
	}

	public String getScoreFile() {
		return scoreFile;
	}

	public String getScoreFileGZ() {
		return scoreFileGZ;
	}

	public void setScoreFileGZ(String sFileGZ) {
		this.scoreFileGZ = sFileGZ;
	}

	public void setGreedy() {
		isGreedy = true;
	}

	public boolean isGreedy() {
		return isGreedy;
	}

	public ExSNPCommandArguments getExSNPCommandArguments() {
		return exSNPCommandArguments;
	}

	public void setExSNPCommandArguments(ExSNPCommandArguments ExSNPCommandArguments) {
		this.exSNPCommandArguments = ExSNPCommandArguments;
	}

	public ProfileCommandArguments getProfileCommandArguments() {
		return profileCommandArguments;
	}

	public void setProfileCommandArguments(ProfileCommandArguments profileCommandArguments) {
		this.profileCommandArguments = profileCommandArguments;
	}

	private String batch;
	private ArrayList<String> bedFiles = NewIt.newArrayList();
	private boolean hasHeader = true;
	private String scoreFile;
	private String scoreFileGZ;
	private boolean isGreedy = false;

	private ExSNPCommandArguments exSNPCommandArguments;
	private ProfileCommandArguments profileCommandArguments;
}
