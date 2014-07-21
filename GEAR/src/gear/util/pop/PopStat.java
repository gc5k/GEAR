package gear.util.pop;

import java.io.IOException;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.BinomialDistribution;
import org.apache.commons.math.distribution.BinomialDistributionImpl;
import org.apache.commons.math.random.RandomDataImpl;

import gear.CmdArgs;
import gear.ConstValues;
import gear.family.popstat.GenotypeMatrix;
import gear.util.Logger;

public class PopStat
{
	public static double[][] calAlleleFrequency(GenotypeMatrix G, int numMarker)
	{
		//[][0]allele freq; [][1]geno freq; [][2] missing rate
		double[][] allelefreq = new double[numMarker][3];
		for (int i = 0; i < G.getGRow(); i++)
		{
			for (int j = 0; j < numMarker; j++)
			{
				int[] c = G.getBiAlleleGenotype(i, j);
				allelefreq[j][c[0]]++;
				allelefreq[j][c[1]]++;
			}
		}

		for (int i = 0; i < numMarker; i++)
		{
			double wa = allelefreq[i][0] + allelefreq[i][1];
			double a = allelefreq[i][0] + allelefreq[i][1] + allelefreq[i][2];
			if (wa > 0)
			{
				for (int j = 0; j < allelefreq[i].length - 1; j++)
				{
					allelefreq[i][j] /= wa;
				}
				allelefreq[i][2] /= a;
			}
			else
			{
				allelefreq[i][0] = Double.NaN;
				allelefreq[i][0] = Double.NaN;
				allelefreq[i][2] = 1;
			}
/*			
			if (allelefreq[i][1] <= 0.5)
			{
				if (allelefreq[i][1] < CmdArgs.INSTANCE.maf_range[0]
						|| allelefreq[i][1] > CmdArgs.INSTANCE.maf_range[1])
				{
					allelefreq[i][1] = 0;
				}
			} 
			else
			{
				if ((1 - allelefreq[i][1]) < CmdArgs.INSTANCE.maf_range[0]
						|| (1 - allelefreq[i][1]) > CmdArgs.INSTANCE.maf_range[1])
				{
					allelefreq[i][1] = 0;
				}
			}
*/
		}

		return allelefreq;
	}

	public static void Imputation (GenotypeMatrix G)
	{
		Logger.printUserLog("Implementing naive imputatation......");
		double[][] f = calAlleleFrequency(G, G.getNumMarker());
		Logger.printUserLog("Missing genotypes will be imputed according to Hardy-Weinberg proportion for each locus with its estimated allele frequency.");

		RandomDataImpl rnd = new RandomDataImpl();
		rnd.reSeed(CmdArgs.INSTANCE.simuSeed);

		int cn = 0;
		for (int i = 0; i < G.getNumMarker(); i++)
		{

			for (int j = 0; j < G.getNumIndivdial(); j++)
			{
				int genoValue = G.getAdditiveScore(j, i);
				if (genoValue == ConstValues.BINARY_MISSING_GENOTYPE)
				{
					cn++;
					int v;
					try
					{
						v = rnd.nextBinomial(2, 1-f[i][0]);
						G.setAdditiveScore(j, i, v);
					}
					catch (MathException e)
					{
						Logger.handleException(e, "Error in imputation.");
					}
				}
			}
		}
		Logger.printUserLog("In total " + cn + " genotypes have been imputed.");
	}
}
