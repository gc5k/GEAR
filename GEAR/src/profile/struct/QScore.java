package profile.struct;

import gear.Parameter;

public class QScore {

	private String delim="\\s+";
	private String SNP;
	private String qscore;
	private boolean isMissing;

	public QScore(String l) {
		String s[] = l.split(delim);
		SNP = s[0];
		qscore = s[1];
		if (Parameter.INSTANCE.isNA(s[1])) {
			isMissing = true;
		} else {
			isMissing = false;
		}
	}

	public String getSNP() {
		return SNP;
	}

	public double getQScore() {
		return Double.parseDouble(qscore);
	}

	public boolean isMissing() {
		return isMissing;
	}

}
