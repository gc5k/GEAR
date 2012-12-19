package profile.struct;

import parameter.Parameter;

public class ScoreUnit {
	
	private String delim="\\s+";
	private String SNP;
	private String RefAllele;
	private String score;
	private boolean isMissing;

	public ScoreUnit(String l) {
		String s[] = l.split(delim);
		SNP = s[0];
		RefAllele = s[1];
		score = s[2];
		if(Parameter.isNA(s[2])) {
			isMissing = true;
		} else {
			isMissing = false;
		}
	}
	
	public ScoreUnit(String snp, String ra, String s) {
		SNP = snp;
		RefAllele = ra;
		score = s;
	}

	public String getSNP() {
		return SNP;
	}
	
	public String getRefAllele() {
		return RefAllele;
	}
	
	public double getScore() {
		return Double.parseDouble(score);
	}
	
	public boolean isMissing() {
		return isMissing;
	}

}
