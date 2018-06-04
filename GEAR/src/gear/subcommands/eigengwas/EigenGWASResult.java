package gear.subcommands.eigengwas;

import java.text.DecimalFormat;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistributionImpl;

import gear.util.Logger;
import gear.util.pop.PopStat;
import gear.util.stat.PrecisePvalue;
import gear.ConstValues;
import gear.family.pedigree.file.SNP;

public class EigenGWASResult {
	private SNP snp;
	private double b;
	private double b_se;
	private double freq;
	private double Z;
	private double P;
	private double PGC;
	private double n1;
	private double freq1;
	private double n2;
	private double freq2;
	private double fst;
	private static DecimalFormat df = new DecimalFormat("0.0000");
	private static DecimalFormat dfE = new DecimalFormat("0.00E000");
	private boolean isGoodLocus;

	private static NormalDistributionImpl unitNormal = new NormalDistributionImpl(0.0D, 1.0D);

	public EigenGWASResult(SNP snp, double freq, double b, double b_se, double n1, double freq1, double n2,
			double freq2, boolean isGood) {
		this.snp = snp;
		this.freq = freq;
		this.b = b;
		this.b_se = b_se;
		Z = b / b_se;
		this.P = getP(Z);

		this.n1 = n1;
		this.freq1 = freq1;
		this.n2 = n2;
		this.freq2 = freq2;
		isGoodLocus = isGood;
		if (isGood) {
			this.fst = PopStat.Fst(freq, (int) n1, (int) n2, freq1, freq2); 
		} else {
			this.fst = Double.NaN;
		}

	}

	public double GetP() {
		return this.P;
	}

	public String printEGWASResult(double gc) {
		StringBuffer sb = new StringBuffer();
		double z1 = Z / Math.sqrt(gc);
		PGC = getP(z1);
		sb.append(snp.getName() + "\t" + snp.getChromosome() + "\t" + snp.getPosition() + "\t" + snp.getFirstAllele()
				+ "\t" + snp.getSecAllele() + "\t" + df.format(freq) + "\t");

		if (!isGoodLocus) {
			sb.append(ConstValues.NA_String + "\t");
		} else {
			if (Math.abs(b) > 0.0001) {
				sb.append(df.format(b) + "\t");
			} else {
				sb.append(dfE.format(b) + "\t");
			}
		}

		if (!isGoodLocus) {
			sb.append(ConstValues.NA_String + "\t");
		} else {
			if (Math.abs(b_se) > 0.0001) {
				sb.append(df.format(b_se) + "\t");
			} else {
				sb.append(dfE.format(b_se) + "\t");
			}
		}

		if (!isGoodLocus) {
			sb.append(ConstValues.NA_String + "\t");
		} else {
			double chi = b * b / (b_se * b_se);

			if (chi > 0.001) {
				sb.append(df.format(chi) + "\t");
			} else {
				sb.append(dfE.format(chi) + "\t");
			}
		}

		if (!isGoodLocus) {
			sb.append(ConstValues.NA_String + "\t");
		} else {
			if (P > 0.0001) {
			sb.append(df.format(P) + "\t");
			} else {
			sb.append(dfE.format(P) + "\t");
			}
		}

		if (!isGoodLocus) {
			sb.append(ConstValues.NA_String + "\t");			
		} else {
			if (PGC > 0.0001) {
				sb.append(df.format(PGC) + "\t");
			} else {
				sb.append(dfE.format(PGC) + "\t");
			}
		}

		sb.append((int) n1 + "\t");

		if (Math.abs(freq1) > 0.0001) {
			sb.append(df.format(freq1) + "\t");
		} else {
			sb.append(dfE.format(freq1) + "\t");
		}

		sb.append((int) n2 + "\t");

		if (Math.abs(freq2) > 0.0001) {
			sb.append(df.format(freq2) + "\t");
		} else {
			sb.append(dfE.format(freq2) + "\t");
		}

		if (!isGoodLocus) {
			sb.append(ConstValues.NA_String + "\t");			
		} else {
			if (Math.abs(fst) > 0.0001) {
				sb.append(df.format(fst) + "\t");
			} else {
			sb.append(dfE.format(fst) + "\t");
			}
		}

		return sb.toString();
	}

	public String printEGWASTab(double gc) {
		StringBuffer sb = new StringBuffer();
		sb.append(snp.getName() + "\t" +  snp.getFirstAllele()
				+ "\t" );

		if (!isGoodLocus) {
			sb.append(ConstValues.NA_String + "\t");
		} else {
			if (Math.abs(b) > 0.0001) {
				sb.append(df.format(b) + "\t");
			} else {
				sb.append(dfE.format(b) + "\t");
			}
		}
		return sb.toString();
	}

	private double getP(double z) {
		double p = 1;
		try {
			if (Math.abs(z) < 8.0D) {
				p = (1.0D - unitNormal.cumulativeProbability(Math.abs(z))) * 2.0D;
			} else {
				p = PrecisePvalue.TwoTailZcumulativeProbability(Math.abs(z));
			}
		} catch (MathException e) {
			Logger.printUserError(e.toString());
		}
		return p;
	}

}
