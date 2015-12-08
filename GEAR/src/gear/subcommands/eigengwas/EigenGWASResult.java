package gear.subcommands.eigengwas;

import java.text.DecimalFormat;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistributionImpl;

import gear.util.Logger;
import gear.util.stat.PrecisePvalue;
import gear.family.pedigree.file.SNP;

public class EigenGWASResult 
{
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
   private static DecimalFormat fmt = new DecimalFormat("0.0000");
   private static DecimalFormat fmtP = new DecimalFormat("0.000E000");

   private static NormalDistributionImpl unitNormal = new NormalDistributionImpl(0.0D, 1.0D);

   public EigenGWASResult(SNP snp, double freq, double b, double b_se, double n1, double freq1, double n2, double freq2, double fst)
   {
	   this.snp = snp;
	   this.freq = freq;
	   this.b = b;
	   this.b_se = b_se;
	   Z = b/b_se;
	   this.P = getP(Z);

	   this.n1 = n1;
	   this.freq1 = freq1;
	   this.n2 = n2;
	   this.freq2 = freq2;
	   this.fst = fst;
   }

   public double GetP()
   {
	   return this.P;
   }

   public String printEGWASResult(double gc)
   {
	   StringBuffer sb = new StringBuffer();
	   double z1 = Z/Math.sqrt(gc);
	   PGC = getP(z1);
       sb.append(snp.getName()+"\t" + snp.getChromosome() + "\t" + snp.getPosition() + "\t" +  snp.getFirstAllele() + "\t" + snp.getSecAllele() + "\t" + freq+ "\t"+fmt.format(b) + "\t" + fmt.format(b_se) + "\t" + fmt.format(Z*Z) + "\t" +fmtP.format(P) + "\t" + fmtP.format(PGC) + "\t"+  (int)n1 + "\t" + fmt.format(freq1) + "\t" + (int)n2 + "\t" + fmt.format(freq2) + "\t" + fmt.format(fst));
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
