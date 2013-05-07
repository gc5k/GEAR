package gear.util.pop;

import family.popstat.GenotypeMatrix;
import gear.CmdArgs;

public class PopStat
{
	public static double[][] calAlleleFrequency(GenotypeMatrix G, int numMarker)
	{
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
		}
		return allelefreq;
	}

}
