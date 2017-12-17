package gear.subcommands.eigengwasepi;

import java.text.DecimalFormat;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistributionImpl;

import gear.family.pedigree.file.SNP;
import gear.util.Logger;
import gear.util.stat.PrecisePvalue;

public class EigenGWASEpiResult
{
	   private SNP snp1;
	   private double Freq1;

	   private double b1;
	   private double b1_se;
	   private double aP1;

	   private double d1;
	   private double d1_se;
	   private double dP1;
	   private double fst1;

	   private SNP snp2;
	   private double Freq2;

	   private double b2;
	   private double b2_se;
	   private double aP2;

	   private double d2;
	   private double d2_se;
	   private double dP2;
	   private double fst2;
	   
	   private double aa;
	   private double aa_se;
	   private double aaP;

	   private static DecimalFormat df = new DecimalFormat("0.0000");
	   private static DecimalFormat dfE = new DecimalFormat("0.00E000");

	   private static NormalDistributionImpl unitNormal = new NormalDistributionImpl(0.0D, 1.0D);

	   public EigenGWASEpiResult(SNP snp1, double Freq1, double fst1, double b1, double b1_se, double d1, double d1_se, 
			   SNP snp2, double Freq2, double fst2, double b2, double b2_se, double d2, double d2_se, double aa, double aa_se)
	   {
		   this.snp1 = snp1;
		   this.Freq1 = Freq1;
		   this.b1 = b1;
		   this.b1_se = b1_se;
		   this.aP1 = getP(b1/b1_se);

		   this.d1 = d1;
		   this.d1_se = d1_se;
		   this.dP1 = getP(d1/d1_se);
		   this.fst1 = fst1;

		   this.snp2 = snp2;
		   this.Freq2 = Freq2;
		   this.b2 = b2;
		   this.b2_se = b2_se;
		   this.aP2 = getP(b2/b2_se);

		   this.d2 = d2;
		   this.d2_se = d2_se;
		   this.dP2 = getP(d2/d2_se);
		   this.fst2 = fst2;

		   this.aa = aa;
		   this.aa_se = aa_se;
		   this.aaP = getP(aa/aa_se);
	   }

	   public double GetAP1()
	   {
		   return this.aP1;
	   }

	   public double GetAP2()
	   {
		   return this.aP2;
	   }

	   public double GetDP1()
	   {
		   return this.dP1;
	   }

	   public double GetDP2()
	   {
		   return this.dP2;
	   }

	   public double GetAAP()
	   {
		   return this.aaP;
	   }

	   public String printEGWASEpiResult()
	   {
		   StringBuffer sb = new StringBuffer();

	       sb.append(snp1.getName()+"\t" + snp1.getChromosome() + "\t" + snp1.getPosition() + "\t" +  snp1.getFirstAllele() + "\t" + snp1.getSecAllele() + "\t" + df.format(Freq1)+ "\t");
	       if (Math.abs(b1) > 0.0001)
	       {
	    	   sb.append(df.format(b1) + "\t");
	       }
	       else
	       {
	    	   sb.append(dfE.format(b1) + "\t");    	   
	       }
	       
	       if (Math.abs(b1_se) > 0.0001)
	       {
	    	   sb.append(df.format(b1_se) + "\t");
	       }
	       else
	       {
	    	   sb.append(dfE.format(b1_se) + "\t");    	   
	       }
	       
	       double chiA1 = b1 * b1 / (b1_se * b1_se);

	       if (chiA1 > 0.001)
	       {
	    	   sb.append(df.format(chiA1) + "\t");
	       }
	       else
	       {
	    	   sb.append(dfE.format(chiA1) + "\t");
	       }

	       if (aP1 > 0.0001)
	       {
	    	   sb.append(df.format(aP1) + "\t");
	       }
	       else
	       {
	    	   sb.append(dfE.format(aP1) + "\t");    	   
	       }

	       if(d1 > 0.0001)
	       {
	    	   sb.append(df.format(d1) + "\t");
	       }
	       else
	       {
	    	   sb.append(dfE.format(d1) + "\t");
	       }

	       if (d1_se > 0.0001)
	       {
	    	   sb.append(df.format(d1_se) + "\t");
	       }
	       else
	       {
	    	   sb.append(dfE.format(d1_se) + "\t");
	       }

	       double chiD1 = d1 * d1 / (d1_se * d1_se);
	       if (chiD1 > 0.0001)
	       {
	    	   sb.append(df.format(chiD1) + "\t");
	       }
	       else
	       {
	    	   sb.append(dfE.format(chiD1) + "\t");
	       }

	       if (dP1 > 0.0001)
	       {
	    	   sb.append(df.format(dP1) + "\t");
	       }
	       else
	       {
	    	   sb.append(dfE.format(dP1) + "\t");
	       }

	       if(Math.abs(fst1) > 0.0001)
	       {
	    	   sb.append(df.format(fst1) + "\t");
	       }
	       else
	       {
	    	   sb.append(dfE.format(fst1) + "\t");    	   
	       }

	       sb.append(snp2.getName()+"\t" + snp2.getChromosome() + "\t" + snp2.getPosition() + "\t" +  snp2.getFirstAllele() + "\t" + snp2.getSecAllele() + "\t" + df.format(Freq2)+ "\t");
	       if(Math.abs(b2) > 0.0001)
	       {
	    	   sb.append(df.format(b2) + "\t");
	       }
	       else
	       {
	    	   sb.append(dfE.format(b2) + "\t");    	   
	       }
	       
	       if(Math.abs(b2_se) > 0.0001)
	       {
	    	   sb.append(df.format(b2_se) + "\t");
	       }
	       else
	       {
	    	   sb.append(dfE.format(b2_se) + "\t");    	   
	       }
	       
	       double chiA2 = b2 * b2 / (b2_se * b2_se);

	       if(chiA2 > 0.001)
	       {
	    	   sb.append(df.format(chiA2) + "\t");
	       }
	       else
	       {
	    	   sb.append(dfE.format(chiA2) + "\t");
	       }

	       if(aP2 > 0.0001)
	       {
	    	   sb.append(df.format(aP2) + "\t");
	       }
	       else
	       {
	    	   sb.append(dfE.format(aP2) + "\t");    	   
	       }

	       if(d2 > 0.0001)
	       {
	    	   sb.append(df.format(d2) + "\t");
	       }
	       else
	       {
	    	   sb.append(dfE.format(d2) + "\t");
	       }

	       if (d2_se > 0.0001)
	       {
	    	   sb.append(df.format(d2_se) + "\t");
	       }
	       else
	       {
	    	   sb.append(dfE.format(d2_se) + "\t");
	       }

	       double chiD2 = d2 * d2/(d2_se * d2_se);
	       
	       if (chiD2 > 0.0001)
	       {
	    	   sb.append(df.format(chiD2) + "\t");
	       }
	       else
	       {
	    	   sb.append(dfE.format(chiD2) + "\t");
	       }

	       if (dP2 > 0.0001)
	       {
	    	   sb.append(df.format(dP2) + "\t");
	       }
	       else
	       {
	    	   sb.append(dfE.format(dP2) + "\t");
	       }

	       if (Math.abs(fst2) > 0.0001)
	       {
	    	   sb.append(df.format(fst2) + "\t");
	       }
	       else
	       {
	    	   sb.append(dfE.format(fst2) + "\t");    	   
	       }

	       if (Math.abs(aa) > 0.0001)
	       {
	    	   sb.append(df.format(aa) + "\t");
	       }
	       else
	       {
	    	   sb.append(dfE.format(aa) + "\t");
	       }

	       if (aa_se > 0.0001)
	       {
	    	   sb.append(df.format(aa_se) + "\t");
	       }
	       else
	       {
	    	   sb.append(dfE.format(aa_se) + "\t");
	       }

	       double chiAA = aa * aa/(aa_se * aa_se);
	       if (chiAA > 0.0001)
	       {
	    	   sb.append(df.format(chiAA) + "\t");
	       }
	       else
	       {
	    	   sb.append(dfE.format(chiAA) + "\t");
	       }

	       if (Math.abs(aaP) > 0.0001)
	       {
	    	   sb.append(df.format(aaP));
	       }
	       else
	       {
	    	   sb.append(dfE.format(aaP));
	       }

	       return sb.toString();
	   }

	   private double getP(double z)
	   {
		   double p = 1;
	       try
	       {
	         if (Math.abs(z) < 8.0D)
	         {
	           p = (1.0D - unitNormal.cumulativeProbability(Math.abs(z))) * 2.0D;
	         }
	         else
	         {
	           p = PrecisePvalue.TwoTailZcumulativeProbability(Math.abs(z));
	         }
	       }
	       catch (MathException e)
	       {
	         Logger.printUserError(e.toString());
	       }
	       return p;
	   }


}
