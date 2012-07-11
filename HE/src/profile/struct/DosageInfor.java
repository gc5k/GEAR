package profile.struct;

import util.SNPMatch;

public class DosageInfor {
	private String delim = "\\s+";
	private String SNP;
	private String RefAllele;
	private String SecAllele;
	private double MAF;
	private double quality;
	private double Rsq;
	private boolean isATGCLocus;
	
	public DosageInfor(String l) {
		String s[] = l.split(delim);
		SNP = s[0];
		RefAllele = s[1];
		SecAllele = s[2];
		isATGCLocus = SNPMatch.Confusion(RefAllele, SecAllele);
		MAF = Double.parseDouble(s[3]);
		quality = Double.parseDouble(s[4]);
		Rsq = Double.parseDouble(s[5]);
	}
	
	public String getSNP() {
		return SNP;
	}
	
	public String getRefAllele() {
		return RefAllele;
	}
	
	public String getSecAllele() {
		return SecAllele;
	}
	
	public double getMAF() {
		return MAF;
	}
	
	public double getQuality() {
		return quality;
	}
	
	public double getRsq() {
		return Rsq;
	}
	
	public boolean isATGCLocus() {
		return isATGCLocus;
	}

}
