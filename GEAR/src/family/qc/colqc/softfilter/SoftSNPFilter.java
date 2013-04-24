package family.qc.colqc.softfilter;

import java.util.HashSet;

import family.popstat.AlleleFrequency;
import family.qc.colqc.SNPFilterInterface;
import gear.CmdArgs;
import gear.util.Logger;
import gear.util.NewIt;

public class SoftSNPFilter
{

	private HashSet<Integer> qualifiedSNPSet;

	private AlleleFrequency al;
	private int[] wseq;
	private int[][] wseq2;
	private int[] bgSNP;

	public SoftSNPFilter(SNPFilterInterface snpFilter, AlleleFrequency al)
	{

		this.al = al;
		this.bgSNP = snpFilter.getBgSeq();

		qualifiedSNPSet = NewIt.newHashSet();
		wseq = snpFilter.getWSeq();
		wseq2 = snpFilter.getWSeq2();

	}

	public void Filter()
	{
		if (wseq != null)
		{
			double[][] allelefreq = al.getAlleleFrequency();
			int count = 0;
			for (int i = 0; i < wseq.length; i++)
			{
				int idx = wseq[i];
				boolean flag = true;
				if (CmdArgs.INSTANCE.genoFlag)
				{
					if (allelefreq[idx][2] > CmdArgs.INSTANCE.geno)
					{
						flag = false;
						continue;
					}
				}
				double f = allelefreq[idx][0] < allelefreq[idx][1] ? allelefreq[idx][0]
						: allelefreq[idx][1];
				if (CmdArgs.INSTANCE.mafFlag)
				{
					if (f < CmdArgs.INSTANCE.maf)
					{
						flag = false;
						continue;
					}
				}
				if (CmdArgs.INSTANCE.maxmafFlag)
				{
					if (f > CmdArgs.INSTANCE.max_maf)
					{
						flag = false;
						continue;
					}
				}
				if (flag)
				{
					qualifiedSNPSet.add(new Integer(idx));
					count++;
				}
			}

			Logger.printUserLog("Number of selected markers: " + count);

			wseq = new int[qualifiedSNPSet.size()];
			int c = 0;
			for (Integer I : qualifiedSNPSet)
			{
				wseq[c++] = I.intValue();
			}
		}
	}

	public int[] getWSeq()
	{
		return wseq;
	}

	public int[][] getWSeq2()
	{
		return wseq2;
	}

	public int[] getBgSeq()
	{
		return bgSNP;
	}

}
