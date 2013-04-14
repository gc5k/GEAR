package gear.profile.struct;

import gear.Parameter;
import gear.util.Logger;

public class ScoreUnit {
	
	private String SNP;
	private String RefAllele;
	private double score;
	private boolean isMissing;
	
	public ScoreUnit(String snp, String ra) {
		SNP = snp;
		RefAllele = ra;
		isMissing = true;
	}
	
	public ScoreUnit(String snp, String ra, double s) {
		SNP = snp;
		RefAllele = ra;
		score = s;
		isMissing = false;
	}

	public String getSNP() {
		return SNP;
	}
	
	public String getRefAllele() {
		return RefAllele;
	}
	
	public double getScore() {
		return score;
	}
	
	public boolean isMissing() {
		return isMissing;
	}

	public static ScoreUnit getNextScoreUnit(gear.util.BufferedReader reader) {
		ScoreUnit scoreUnit = null;
		
		while (true) {
			String tokens[] = reader.readTokens();
			
			if (tokens == null) {
				break;
			}
			
			if (tokens.length != 3) {
				Logger.printUserError(
					"The format of the score file '" + reader.getFileName() + "' is incorrect: " +
					"A score file should consists of 3 columns, but line " + reader.getCurLineNum() + " contains " + tokens.length + " column(s)."
				);
				System.exit(1);
			}
			
			if (Parameter.INSTANCE.isNA(tokens[2])) {
				scoreUnit = new ScoreUnit(/* SNP = */ tokens[0], /* RefAllele = */ tokens[1]);
			}
			else {
				try {
					scoreUnit = new ScoreUnit(/* SNP = */ tokens[0], /* RefAllele = */ tokens[1], /* score = */ Double.parseDouble(tokens[2]));
				}
				catch (NumberFormatException e) {
					Logger.handleException(e,
						"Line " + reader.getCurLineNum() + " of the score file '" + reader.getFileName() + "' " +
					    "contains an error: '" + tokens[2] + "' is not a valid score."
					);
				}
			}
			
			break;
		}
		
		return scoreUnit;
	}

}
