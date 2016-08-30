package gear.subcommands.oath.nss;

import java.text.DecimalFormat;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistributionImpl;

import gear.util.Logger;
import gear.util.stat.PrecisePvalue;
import gear.family.pedigree.file.SNP;

public class NSSGWASResult 
{
   private SNP snp;
   private double b;
   private double b_se;
   private double freq;
   private double vg;
   private double Z;
   private double P;
   private static DecimalFormat df = new DecimalFormat("0.0000");
   private static DecimalFormat dfE = new DecimalFormat("0.00E000");

   private static NormalDistributionImpl unitNormal = new NormalDistributionImpl(0.0D, 1.0D);

   public NSSGWASResult(SNP snp, double freq, double vg, double b, double b_se)
   {
	   this.snp = snp;
	   this.freq = freq;
	   this.vg = vg;
	   this.b = b;
	   this.b_se = b_se;
	   Z = b/b_se;
	   this.P = getP(Z);
   }

   public double GetP()
   {
	   return this.P;
   }

   public String printEGWASResult(double gc)
   {
	   StringBuffer sb = new StringBuffer();
       sb.append(snp.getName()+"\t" + snp.getChromosome() + "\t" + snp.getPosition() + "\t" +  snp.getFirstAllele() + "\t" + snp.getSecAllele() + "\t" + df.format(freq)+ "\t" + df.format(vg) + "\t");
       if(Math.abs(b) > 0.0001)
       {
    	   sb.append(df.format(b) + "\t");
       }
       else
       {
    	   sb.append(dfE.format(b) + "\t");    	   
       }
       
       if(Math.abs(b_se) > 0.0001)
       {
    	   sb.append(df.format(b_se) + "\t");
       }
       else
       {
    	   sb.append(dfE.format(b_se) + "\t");    	   
       }
       
       double chi = b * b / (b_se * b_se);

       if(chi > 0.001)
       {
    	   sb.append(df.format(chi) + "\t");
       }
       else
       {
    	   sb.append(dfE.format(chi) + "\t");
       }

       if(P > 0.0001)
       {
    	   sb.append(df.format(P) + "\t");
       }
       else
       {
    	   sb.append(dfE.format(P) + "\t");    	   
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
