package gear.subcommands.mlmmeta;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.Set;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;
import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.EigenDecompositionImpl;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;

import gear.gwassummary.GWASReader;
import gear.gwassummary.MetaStat;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.subcommands.mlmmeta.util.MLMCovMatrix;
import gear.subcommands.mlmmeta.util.MLMRes;
import gear.util.BufferedReader;
import gear.util.Logger;
import gear.util.NewIt;
import gear.util.SNPMatch;
import gear.util.stat.MLM;

public class MLMMetaCommandImpl extends CommandImpl
{

	@Override
	public void execute(CommandArguments cmdArgs)
	{
		mlmMetaArgs = (MLMMetaCommandArguments) cmdArgs;

		if (mlmMetaArgs.isQT())
		{
			Logger.printUserLog("Analysing summary statistics analysis for quantitative traits.\n");			
		}
		else
		{
			Logger.printUserLog("Analysing summary statistics analysis for case-contrl studies.\n");			
		}

		FileKeep = new boolean[mlmMetaArgs.getMetaFile().length];
		Arrays.fill(FileKeep, true);

		if (mlmMetaArgs.IsKeepFile() || mlmMetaArgs.IsRevFile())
		{
			FilterFiles();
		}

		if (mlmMetaArgs.isGC())
		{
			Logger.printUserLog("Genomic-control only applies for cohorts that have gc factor > 1.");
		}
		
		if (mlmMetaArgs.isMLM())
		{
			generateMLMMatrix();			
		}
		generateCorMatrix();

		gReader = new GWASReader(mlmMetaArgs.getMetaFile(), FileKeep, mlmMetaArgs.getKeys(), mlmMetaArgs.isQT(), mlmMetaArgs.isGZ(), 
				mlmMetaArgs.getChr() != null? true:false,
			    mlmMetaArgs.getChr() != null? mlmMetaArgs.getChr():mlmMetaArgs.getNotChr());
		gReader.Start(false);

		if (gReader.getNumMetaFile() < 2)
		{
			Logger.printUserError("At least two summary statistic files should be specified.\n");
			Logger.printUserError("GEAR quitted.\n");
			System.exit(0);
		}

		MetaAnalysis();
	}

	private void generateMLMMatrix()
	{
		ArrayList<String> tWorkingMetaFile = NewIt.newArrayList();

		int cn = 0;
		for(int i = 0; i < FileKeep.length; i++)
		{
			if(FileKeep[i]) cn++;
		}
		mlmMat = new double[cn][cn];

		if (!mlmMetaArgs.isMLM())
		{
			Logger.printUserLog("No mlm matrix is specified. The default correlation (digonal matrix) will be used.");
		}
		else
		{
			int NumMetaFile = mlmMetaArgs.getMetaFile().length;
			BufferedReader bf = BufferedReader.openTextFile(mlmMetaArgs.getMLMFile(), "mlm file.");
			Logger.printUserLog("Reading '" + mlmMetaArgs.getMLMFile() + "'.");

			String[] d = null;
			int cnt = 0;
			int cIdx = 0;
			while ( (d = bf.readTokens())!= null )
			{
				if (d.length != NumMetaFile )
				{
					Logger.printUserError("incorrect '" + mlmMetaArgs.getMLMFile() + "'.");
					System.exit(0);
				}
				if(FileKeep[cnt])
				{
					int c = 0;
					for(int i = 0; i < d.length; i++)
					{
						if(FileKeep[i])
						{
							mlmMat[cIdx][c++] = Double.parseDouble(d[i]);
						}
					}
					cIdx++;
				}
				cnt++;
			}

			for(int i = 0; i < mlmMat.length; i++)
			{
				for(int j = 0; j < i; j++)
				{
					mlmMat[j][i] = mlmMat[i][j];
				}
			}

			Logger.printUserLog(mlmMat.length + "X" + mlmMat.length + " mlm matrix has been read in.");

			for(int i = 0; i < FileKeep.length; i++)
			{
				if(FileKeep[i])
				{
					tWorkingMetaFile.add(mlmMetaArgs.getMetaFile()[i]);
				}
			}
		}
	}

	private void generateCorMatrix()
	{
		ArrayList<String> tWorkingMetaFile = NewIt.newArrayList();

		int cn = 0;
		for(int i = 0; i < FileKeep.length; i++)
		{
			if(FileKeep[i]) cn++;
		}
		corMat = new double[cn][cn];
		zMat = new double[cn][cn];

		if (mlmMetaArgs.getCMFile() == null)
		{
			for(int i = 0; i < corMat.length; i++)
			{
				corMat[i][i] = 1;
			}
			Logger.printUserLog("No correlation matrix is specified. The default correlation (digonal matrix) will be used.");
		}
		else
		{
			int NumMetaFile = mlmMetaArgs.getMetaFile().length;
			BufferedReader bf = BufferedReader.openTextFile(mlmMetaArgs.getCMFile(), "cm file.");
			Logger.printUserLog("Reading '" + mlmMetaArgs.getCMFile() + "'.");

			String[] d = null;
			int cnt = 0;
			int cIdx = 0;
			while ( (d = bf.readTokens())!= null )
			{
				if (d.length != NumMetaFile )
				{
					Logger.printUserError("incorrect '" + mlmMetaArgs.getCMFile() + "'.");
					System.exit(0);
				}
				if(FileKeep[cnt])
				{
					int c = 0;
					for(int i = 0; i < d.length; i++)
					{
						if(FileKeep[i])
						{
							corMat[cIdx][c++] = Double.parseDouble(d[i]);
						}
					}
					cIdx++;
				}
				cnt++;
			}

			for(int i = 0; i < corMat.length; i++)
			{
				for(int j = 0; j < i; j++)
				{
					zMat[i][j] = corMat[j][i];
					zMat[j][i] = corMat[j][i];
					corMat[j][i] = corMat[i][j];
				}
			}

			Logger.printUserLog(corMat.length + "X" + corMat.length + " correlation matrix has been read in.");
			
			for(int i = 0; i < FileKeep.length; i++)
			{
				if(FileKeep[i])
				{
					tWorkingMetaFile.add(mlmMetaArgs.getMetaFile()[i]);
				}
			}

			if (mlmMetaArgs.getDiag())
			{
				RemMetaIdx = Zprune(tWorkingMetaFile);
				if(RemMetaIdx != null)
				{
					for(int i = 0; i < RemMetaIdx.length; i++)
					{
						FileKeep[RemMetaIdx[i]] = false;
					}		
				}
			}
			if(mlmMetaArgs.getNaive())
			{
				Logger.printUserLog("Force the " + corMat.length + "X" + corMat.length + "correlation matrix to be diagonal matrix for naive meta-analysis.");
				corMat = new double[corMat.length][corMat.length];
				for(int i = 0; i < corMat.length; i++)
				{
					corMat[i][i] = 1;
				}
			}
		}
	}

	private int[] Zprune(ArrayList<String> workingMetaFile)
	{
		Logger.printUserLog("Starting matrix pruning...");
		NormalDistribution nd = new NormalDistributionImpl();
		double threshold = 0;
		try
		{
			threshold = nd.inverseCumulativeProbability(1-0.05/(zMat.length * (zMat.length - 1) / 2));
		}
		catch (MathException e)
		{
			e.printStackTrace();
		}
		
		ArrayList<Integer> Zidx = NewIt.newArrayList();
		ArrayList<Integer> RemIdx = NewIt.newArrayList();
		for(int i = 0; i < zMat.length; i++)
		{
			Zidx.add(i);
		}

		double[][] zm = new double[Zidx.size()][Zidx.size()];
		double[][] cm = new double[Zidx.size()][Zidx.size()];
		for(int i = 0; i < Zidx.size(); i++)
		{
			for(int j = 0; j < Zidx.size(); j++)
			{
				zm[i][j] = zMat[Zidx.get(i).intValue()][Zidx.get(j).intValue()];
				cm[i][j] = corMat[Zidx.get(i).intValue()][Zidx.get(j).intValue()];
			}
		}

		RealMatrix gg = new Array2DRowRealMatrix(cm);

		boolean isNonSingular = (new LUDecompositionImpl(gg)).getSolver().isNonSingular();
		System.out.println(isNonSingular);
		EigenDecompositionImpl EI= new EigenDecompositionImpl(gg, 0.00000001);
		double[] ei = EI.getRealEigenvalues();
		Arrays.sort(ei);
		int it = 0;
		while(ei[0] < 0)
		{
			int max = 0;
			int maxIdx = 0;
			int[] cnt = new int[zm.length];
			for (int i = 0; i < zm.length; i++)
			{
				for (int j = 0; j < zm[i].length; j++)
				{
					if (threshold < zm[i][j])
					{
						cnt[i]++;
					}
				}
				if (max < cnt[i])
				{
					max = cnt[i];
					maxIdx = i;
				}
			}

			Logger.printUserLog("Remove " + workingMetaFile.get(Zidx.get(maxIdx)));
			RemIdx.add(Zidx.get(maxIdx));
			Zidx.remove(maxIdx);

			zm = new double[Zidx.size()][Zidx.size()];
			cm = new double[Zidx.size()][Zidx.size()];
			for(int i = 0; i < Zidx.size(); i++)
			{
				for(int j = 0; j < Zidx.size(); j++)
				{
					zm[i][j] = zMat[Zidx.get(i).intValue()][Zidx.get(j).intValue()];
					cm[i][j] = corMat[Zidx.get(i).intValue()][Zidx.get(j).intValue()];
				}
			}

			gg = new Array2DRowRealMatrix(cm);

			isNonSingular = (new LUDecompositionImpl(gg)).getSolver().isNonSingular();
			EI= new EigenDecompositionImpl(gg, 0.00000001);
			ei = EI.getRealEigenvalues();
			Arrays.sort(ei);
			it++;
		}

		if (it > 0)
		{
			corMat = new double[cm.length][cm.length];
			for(int i = 0; i < cm.length; i++)
			{
				System.arraycopy(cm[i], 0, corMat[i], 0, cm.length);
			}

			Collections.sort(RemIdx);
			RemMetaIdx = new int[RemIdx.size()];
			for(int i = 0; i < RemIdx.size(); i++)
			{
				System.out.println(RemIdx.get(i).intValue());
				RemMetaIdx[i] = RemIdx.get(i).intValue();
			}
			if (it == 1)
			{
				Logger.printUserLog("Removed " + RemMetaIdx.length + " cohort.");				
			}
			else
			{
				Logger.printUserLog("Removed " + RemMetaIdx.length + " cohorts.");
			}
		}
		else
		{
			Logger.printUserLog("No cohorts have been removed in matrix diagnosis.");
		}

		return RemMetaIdx;
	}

	private void MetaAnalysis()
	{
		Logger.printUserLog("Starting meta-analysis...");
		int singularCnt = 0;
		int atgcCnt = 0;
		Set<String> snps = gReader.getMetaSNPTable().keySet();
		for (Iterator<String> e=snps.iterator(); e.hasNext();)
		{
			String snp = e.next();
			ArrayList<Integer> Int = gReader.getMetaSNPTable().get(snp);
			for(int i = 0; i < (Int.size() - 1); i++)
			{
				if(Int.get(i).intValue() != 0) break; 
			}

			if (mlmMetaArgs.isFullSNPOnly())
			{
				if (Int.get(Int.size()-1).intValue() != (Int.size() -1))
				{
					continue;
				}
			}
			MLMCovMatrix covMat = new MLMCovMatrix(snp, Int, mlmMetaArgs.getQTsize(), mlmMat, corMat, gReader, mlmMetaArgs.isGC(), mlmMetaArgs.isGCALL(), mlmMetaArgs.IsAdjOverlappingOnly());

			MLMRes mlmRes = getBeta(covMat);

			MLM mlm = new MLM(covMat.getCovMatrix(), covMat.getA(), mlmRes.getRawBeta(), true);
			mlm.setSilent(true);
			mlm.MINQUE();
			RealMatrix bt = mlm.getBeta();
			RealMatrix bt_v = mlm.getVarBeta();
			RealMatrix Vc =mlm.getVC();

			mlmRes.SetB(bt.getEntry(0, 0));
			mlmRes.SetSE(Math.sqrt(bt_v.getEntry(0, 0)));
			mlmRes.SetVC(Vc.getEntry(0, 0));
			grArray.add(mlmRes);
		}
		Collections.sort(grArray);

		Logger.printUserLog("In total "+ snps.size() + " loci have been read.");
		Logger.printUserLog("In total "+ grArray.size() + " loci have been used for meta-analysis.");
		if (singularCnt > 0)
		{
			Logger.printUserLog(singularCnt + " loci were excluded from analyais because of singular matrix.");
		}

		if (!mlmMetaArgs.isKeepATGC())
		{
			Logger.printUserLog(atgcCnt + " ambiguous loci have been eliminated.");
		}

		PrintMLMresults();
	}

	private MLMRes getBeta(MLMCovMatrix mlmcovMat)
	{
		String SNP = mlmcovMat.getSNP();
		int[] idx = mlmcovMat.getCohortIdx();
		int cohort = idx.length;

		StringBuffer direction = new StringBuffer();

		for(int i = 0; i < gReader.getCohortNum(); i++)
		{
			direction.append('?');
		}

		double[] beta = new double[idx.length];
		MLMRes gr = new MLMRes(cohort);
		MetaStat ms = null;
		boolean isAmbiguousLocus = false;
		for (int i = 0; i < idx.length; i++)
		{
			char sign = '+';
			float b = 0;
			boolean match = true;
			ms = gReader.getMetaStat().get(idx[i]).get(SNP);
			if (i == 0)
			{
				b = ms.getEffect();
				gr.SetSNP(ms.getSNP());
				gr.SetChr(ms.getChr());
				gr.SetBP(ms.getBP());
				gr.SetA1(ms.getA1());
				gr.SetA2(ms.getA2());
				isAmbiguousLocus = SNPMatch.isAmbiguous(ms.getA1(), ms.getA2());
			}
			else
			{
				match = SNPMatch.isAllelesMatchForTwoLoci(gr.GetA1(), gr.GetA2(), ms.getA1(), ms.getA2());
				if (!match)
				{
					match = SNPMatch.isAllelesFlipMatchForTwoLoci(gr.GetA1(), gr.GetA2(), ms.getA1(), ms.getA2());
				}
				
				if (match)
				{
					b = ms.getEffect();

					if (gr.GetA1() == ms.getA1() || gr.GetA1() == SNPMatch.Flip(ms.getA1())) //match A1 in the second meta
					{
						b = 1 * b;
					}
					else if (gr.GetA1() == ms.getA2() || gr.GetA1() == SNPMatch.Flip(ms.getA2())) //match A2 in the second meta
					{
						b = -1 * b;
					}

				}
				else
				{
					sign = ',';
				}
			}

			if (match)
			{
				if(b == 0)
				{
					sign = '0';
				}
				else if(b > 0)
				{
					sign = '+';
				}
				else
				{
					sign = '-';
				}
			}
			direction.setCharAt(idx[i], sign);
			beta[i] = b;
		}
		
		gr.SetAmbi(isAmbiguousLocus);
		gr.SetDirect(direction.toString());
		gr.SetRawBeta(beta);
		return gr;
	}

	private void PrintMLMresults()
	{
        PrintWriter writer = null;
        try
        {
        	writer = new PrintWriter(new BufferedWriter(new FileWriter(mlmMetaArgs.getOutRoot()+".mlm-meta")));
        	Logger.printUserLog("Writting detailed test statistics into '"+mlmMetaArgs.getOutRoot() + ".mlm-meta.'\n");
        }
		catch (IOException e)
		{
			Logger.handleException(e, "An I/O exception occurred when writing '" + mlmMetaArgs.getOutRoot() + ".mlm-meta" + "'.\n");
		}

		for(int i = 0; i < grArray.size(); i++)
		{
			MLMRes gr = grArray.get(i);
			if(i == 0)
			{
				writer.write(gr.printTitle() + "\n");
			}
			writer.write(gr.toString()+ "\n");
		}
		writer.close();
	}

	private void FilterFiles()
	{
		String[] metaF = mlmMetaArgs.getMetaFile();

		if(mlmMetaArgs.IsKeepFile())
		{
			Arrays.fill(FileKeep, false);
			String[] kf = mlmMetaArgs.getKeepCohortFile();
			for(int i = 0; i < metaF.length; i++)
			{
				for(int j = 0; j < kf.length; j++ )
				{
					if(metaF[i].compareTo(kf[j]) == 0)
					{
						FileKeep[i] = true;
					}
				}
			}
		}

		if(mlmMetaArgs.IsRevFile())
		{
			Arrays.fill(FileKeep, true);
			String[] kf = mlmMetaArgs.getRemoveCohortFile();
			for(int i = 0; i < metaF.length; i++)
			{
				for(int j = 0; j < kf.length; j++ )
				{
					if(metaF[i].compareTo(kf[j]) == 0)
					{
						FileKeep[i] = false;
					}
				}
			}
		}

		for(int i = 0; i < FileKeep.length; i++)
		{
			System.out.println(metaF[i] + " " +FileKeep[i]);
		}
	}

	private MLMMetaCommandArguments mlmMetaArgs;
	private GWASReader gReader;
	private double[][] mlmMat;
	private double[][] corMat;
	private double[][] zMat;
	private boolean[] FileKeep;
	private int[] RemMetaIdx;
	
	private ArrayList<MLMRes> grArray = NewIt.newArrayList();
}
