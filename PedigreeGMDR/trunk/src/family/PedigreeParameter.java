package family;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.regex.Pattern;

public class PedigreeParameter {

	protected ArrayList<String> lines;
	protected String ConfigFileName;
	protected BufferedReader buffer;
	protected boolean isLouAlgorithm; // 0
	protected boolean usingFounderGenotype; //1
	protected boolean usingChildrenGenotype; //2
	protected String ped_file; // 3
	protected String phe_file; // 4
	protected String converted_ped_file; // 5
	protected String converted_phe_file; // 6
	protected String id_file; // 7
	protected int[] cov_idx; // 8 starts from 1
	protected int phe_idx; // 9 starts from 1
	protected int scoreBuildMethod; // 10; 1 for linear, 2 for logistic
	protected boolean scoreAdjustmentScheme; // 11
	protected boolean scoreBuildWithFounder; // 12
	protected boolean scoreBuildWithChildren; // 13
	protected int replication; // 14
	protected long seed; // 15

	public PedigreeParameter() {
		lines = new ArrayList();
	}

	public void read(String file) throws IOException {
		ConfigFileName = file;
		buffer = new BufferedReader(new FileReader(new File(file)));
		sweepComments();
		parseValue();
	}

	public void sweepComments() throws IOException {
		boolean flag = true;
		String line;
		while ((line = buffer.readLine()) != null) {
			if (Pattern.matches("^//.*", line)) {// empty line
				continue;
			} else {
				lines.add(line);
			}
		}
	}

	public void parseValue() {
		isLouAlgorithm = Boolean.parseBoolean(lines.get(0));
		usingFounderGenotype = Boolean.parseBoolean(lines.get(1));
		usingChildrenGenotype = Boolean.parseBoolean(lines.get(2));
		ped_file = lines.get(3);
		phe_file = lines.get(4);
		converted_ped_file = lines.get(5);
		converted_phe_file = lines.get(6);
		id_file = lines.get(7);
		String[] c = lines.get(8).split(",");
		cov_idx = new int[c.length];
		for (int i = 0; i < c.length; i++) {
			cov_idx[i] = Integer.parseInt(c[i]) - 1;
		}
		phe_idx = Integer.parseInt(lines.get(9)) - 1;
		scoreBuildMethod = Integer.parseInt(lines.get(10));
		scoreAdjustmentScheme = Boolean.parseBoolean(lines.get(11));
		scoreBuildWithFounder = Boolean.parseBoolean(lines.get(12));
		scoreBuildWithChildren = Boolean.parseBoolean(lines.get(13));
		replication = Integer.parseInt(lines.get(14));
		seed = Long.parseLong(lines.get(15));
	}

	public void setIsLouAlgorithm(boolean Lou) {
		isLouAlgorithm = Lou;
	}
	
	public boolean isLouAlgorithm() {
		return isLouAlgorithm;
	}

	public void setUsingFounderGenotype(boolean Founder) {
		usingFounderGenotype = Founder;
	}
	
	public boolean usingFounderGenotype() {
		return usingFounderGenotype;
	}

	public void setUsingChildrenGenotype(boolean Founder) {
		usingChildrenGenotype = Founder;
	}

	public boolean usingChildrenGenotype() {
		return usingChildrenGenotype;
	}

	public void setPedigreeFile(String ped) {
		ped_file = ped;
	}

	public String getPedigreeFile() {
		return ped_file;
	}
	
	public void setPhenotypeFile(String phe) {
		phe_file = phe;
	}
	
	public String getPhenotypeFile() {
		return phe_file;
	}

	public void setConvertedPedigreeFile(String cPed) {
		converted_ped_file = cPed;
	}

	public String getConvertedPedigreeFile() {
		return converted_ped_file;
	}

	public void setConvertedPhenotypeFile(String cPhe) {
		converted_phe_file = cPhe;
	}
	
	public String getConvertedPhenotypeFile() {
		return converted_phe_file;
	}
	
	public void setFamilyIDFile(String fID) {
		id_file = fID;
	}
	
	public String getFamilyIDFile() {
		return id_file;
	}

	public void setCovariateIndex(String ind) {
		String[] c = ind.split(",");
		cov_idx = new int[c.length];
		for (int i = 0; i < c.length; i++) {
			cov_idx[i] = Integer.parseInt(c[i]) - 1;
		}
	}

	public int[] getCovarianteIndex() {
		return cov_idx;
	}

	public void setPhenotypeIndex(String ind) {
		phe_idx = Integer.parseInt(ind) - 1;
	}

	public int getPhenotypeIndex() {
		return phe_idx;
	}

	public void setScoreBuildMethod(int m) {
		//1 for linear, 2 for logistic
		scoreBuildMethod = m;
	}

	public int getScoreBuildMethod() {
		return scoreBuildMethod;
	}

	public void setAdjustScore(boolean b) {
		scoreAdjustmentScheme = b;
	}

	public boolean isAdjustScore() {
		return scoreAdjustmentScheme;
	}

	public void setScoreBuildWithFounder(boolean b) {
		scoreBuildWithFounder = b;
	}

	public boolean getScoreBuildWithFounder() {
		return scoreBuildWithFounder;
	}

	public void setScoreBuildWithChildren(boolean b) {
		scoreBuildWithChildren = b;
	}

	public boolean getScoreBuildWithChildren() {
		return scoreBuildWithChildren;
	}	
	
	public void setReplication(int r) {
		replication = r;
	}
	
	public int getReplication() {
		return replication;
	}
	
	public void setSeed(long s) {
		seed = s;
	}
	
	public long getSeed() {
		return seed;
	}	
}
