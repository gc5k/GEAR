package family;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.regex.Pattern;

public class PedigreeParameter {

	protected ArrayList<String> lines;
	protected BufferedReader buffer;
	protected String filename;
	protected String ped_file; // 2
	protected String phe_file; // 3
	protected String converted_ped_file; // 4
	protected String converted_phe_file; // 5
	protected String id_file; // 6
	protected int[] cov_idx; // 7 starts from 1
	protected int[] phe_idx; // 8 starts from 1
	protected int scoreBuildMethod; // 9; 1 for linear, 2 for logistic
	protected boolean scoreAdjustmentScheme; // 10
	protected boolean scoreBuildWithFounder; // 11
	protected int replication; // 12

	public PedigreeParameter() {
		lines = new ArrayList();
	}

	public void read(String file) throws IOException {
		filename = file;
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
		ped_file = lines.get(0);
		phe_file = lines.get(1);
		converted_ped_file = lines.get(2);
		converted_phe_file = lines.get(3);
		id_file = lines.get(4);
		String[] c = lines.get(5).split(",");
		cov_idx = new int[c.length];
		for (int i = 0; i < c.length; i++) {
			cov_idx[i] = Integer.parseInt(c[i]) - 1;
		}
		String[] p = lines.get(6).split(",");
		phe_idx = new int[p.length];
		for (int i = 0; i < p.length; i++) {
			phe_idx[i] = Integer.parseInt(p[i]) - 1;
		}
		scoreBuildMethod = Integer.parseInt(lines.get(7));
		scoreAdjustmentScheme = Boolean.parseBoolean(lines.get(8));
		scoreBuildWithFounder = Boolean.parseBoolean(lines.get(9));
		replication = Integer.parseInt(lines.get(10));
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
		String[] c = ind.split(",");
		phe_idx = new int[c.length];
		for (int i = 0; i < c.length; i++) {
			phe_idx[i] = Integer.parseInt(c[i]) - 1;
		}
	}

	public int[] getPhenotypeIndex() {
		return phe_idx;
	}

	public void setScoreBuildMethod(int m) {
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

	public void setReplication(int r) {
		replication = r;
	}
	
	public int getReplication() {
		return replication;
	}
}
