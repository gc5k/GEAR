package gear.subcommands.arch;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;

import gear.ConstValues;
import gear.family.pedigree.file.MapFile;
import gear.family.pedigree.file.SNP;
import gear.family.plink.PLINKParser;
import gear.family.popstat.GenotypeMatrix;
import gear.family.qc.rowqc.SampleFilter;
import gear.gwassummary.GWASReader;
import gear.gwassummary.MetaStat;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.sumstat.qc.rowqc.SumStatQC;
import gear.util.BufferedReader;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;
import gear.util.SNPMatch;
import gear.util.stat.JModel;

public class ArchCommandImpl extends CommandImpl
{

	@Override
	public void execute(CommandArguments cmdArgs) 
	{
		archArgs = (ArchCommandArguments) cmdArgs;
		
		PLINKParser pp = PLINKParser.parse(this.archArgs);

		this.sf = new SampleFilter(pp.getPedigreeData(), pp.getMapData());
		this.ssQC = new SumStatQC(pp.getPedigreeData(), pp.getMapData(), this.sf);
		this.mapFile = this.ssQC.getMapFile();
		this.gm = new GenotypeMatrix(this.ssQC.getSample());
		
		readExtract();
		readSNPinfo();
//		readLDscore();
		lineup();
		jmapping();
		//he();
	}

	private void jmapping()
	{
		for (int i = 0; i < snpExtract.size(); i++)
		{
			int index = snpListAdjName.indexOf(snpExtract.get(i));
			if (index == -1) continue;
			
			SNP snpLead = snpListAdj.get(index);
			
			long start = (long) (snpLead.getPosition() - archArgs.getWindow() * 1000000);
			long end = (long) (snpLead.getPosition() + archArgs.getWindow() * 1000000);

			ArrayList<Integer> snpIdx = NewIt.newArrayList();
			snpIdx.add(index);
			int id1 = index - 1;
			while ( id1 >=0 )
			{
				SNP snp1 = snpListAdj.get(id1);
				if ( (snp1.getChromosome().compareTo(snpLead.getChromosome())) == 0 && snp1.getPosition() > start)
				{
					snpIdx.add(id1);
					id1--;
				}
				else
				{
					break;
				}
			}

			int id2 = index + 1;
			while ( id2 < snpListAdj.size() )
			{
				SNP snp2 = snpListAdj.get(id2);
				if ((snp2.getChromosome().compareTo(snpLead.getChromosome())) == 0 && snp2.getPosition() < end)
				{
					snpIdx.add(id2);
					id2++;
				}
				else
				{
					break;
				}
			}

			Collections.sort(snpIdx);

			Logger.printUserLog(snpIdx.size() + " SNPs found around " + snpLead.getName());
			ArrayList<Integer> ldIdx = NewIt.newArrayList();
			for (int l1 = 0; l1 < snpIdx.size(); l1++)
			{
				System.out.println(snpListAdj.get(snpIdx.get(l1)).getName());
				ldIdx.add(snpSel.get(snpIdx.get(l1)));
			}
			double[][] XX = getLDMat(ldIdx);

			for (int l1 = 0; l1 < XX.length; l1++)
			{
				for (int l2 = 0; l2 < XX.length; l2++)
				{
					System.out.print(XX[l1][l2] + " ");
				}
				System.out.println();
			}
			double[] b = new double[snpIdx.size()];

			for (int k = 0; k < b.length; k++)
			{
				b[k] =  effectSel.get(snpIdx.get(k));
			}

			JModel jm = new JModel(XX, b);
			double[] jb = jm.getJB();
			if (jm.isSingular())
			{
				Logger.printUserLog(snpLead.getName() + " is singular.");
			}
			else
			{
				PrintStream jFile = FileUtil.CreatePrintStream(snpLead.getName() + ".txt");
				for (int k = 0; k < snpIdx.size(); k++)
				{
					jFile.println(snpListAdj.get(snpIdx.get(k)).getName() + " " + jb[k] + " " + effectSel.get(snpIdx.get(k)));
				}
				jFile.close();
			}
		}
	}

	private double[][] getLDMat(ArrayList<Integer> snpIdx)
	{
		double[][] ldMat = new double[snpIdx.size()][snpIdx.size()];
		for (int i = 0; i < snpIdx.size(); i++)
		{
			for (int j = 0; j <= i; j++)
			{

				int cnt = 0;
				double sg1 = 0;
				double sg2 = 0;
				double ssg1 = 0;
				double ssg2 = 0;
				double crs = 0;
				for (int k = 0; k < gm.getNumIndivdial(); k++)
				{
					int g1 = gm.getAdditiveScoreOnFirstAllele(k, snpIdx.get(i));
					int g2 = gm.getAdditiveScoreOnFirstAllele(k, snpIdx.get(j));
					if (g1 == ConstValues.BINARY_MISSING_GENOTYPE || g2 == ConstValues.BINARY_MISSING_GENOTYPE)
					{
						continue;
					}
					sg1 += g1;
					sg2 += g2;
					ssg1 += g1*g1;
					ssg2 += g2*g2;
					crs += g1*g2;
					cnt++;
				}
				if (cnt > 10)
				{
					double m1 = sg1/(2*cnt);
					double m2 = sg2/(2*cnt);
					double v1 = ssg1/(4*cnt) - m1 * m1;
					double v2 = ssg2/(4*cnt) - m2 * m2;
					double cv = crs/(4*cnt) - m1 * m2;
					ldMat[i][j] = ldMat[j][i] = cv/Math.sqrt(v1*v2);
				}
			}
		}
		return ldMat;
	}

//	private void he()
//	{
//		float d1 = 0;
//		float d2 = 0;
//
//		for(int m = 0; m < snpSel.size(); m++)
//		{
//			SNP snp0 = snpList.get(m);
//			for(int l1 = 0; l1 < snpSel.size(); l1++)
//			{
//				SNP snp1 = snpList.get(snpSel.get(l1));
//
//				if (archArgs.isWindow())
//				{
//					if (snp0.getChromosome().compareTo(snp1.getChromosome()) != 0) continue;
//					if (Math.abs(snp0.getPosition() - snp1.getPosition()) > (archArgs.getWindow() * 1000000) ) continue;						
//				}
//
//				for(int l2 = 0; l2 < snpSel.size(); l2++)
//				{
//					SNP snp2 = snpList.get(snpSel.get(l2));
//					
//					if (archArgs.isWindow())
//					{
//						if (snp1.getChromosome().compareTo(snp2.getChromosome()) != 0) continue;
//						if (Math.abs(snp1.getPosition() - snp2.getPosition()) > (archArgs.getWindow() * 1000000) ) continue;						
//					}
//
//					double v = ldMatOrg[snpSel.get(m)][snpSel.get(l1)] * ldMatOrg[snpSel.get(m)][snpSel.get(l2)] * effect.get(l1) * effect.get(l2);
//
//					if (l1 == l2)
//					{
//						d1 += v;
//					}
//					else
//					{
//						d2 += v;
//					}
//				}
//			}
//		}
//
//		d1 /= snpSel.size();
//		d2 /= snpSel.size();
//		
//		Logger.printUserLog("d1 is: " + d1);
//		Logger.printUserLog("d2 is: " + d2);
//
//		double ldSum = 0;
//		for(int l1 = 0; l1 < snpSel.size(); l1++)
//		{
//			for(int l2 = 0; l2 < snpSel.size(); l2++)
//			{
//				ldSum += ldMatOrg[snpSel.get(l1)][snpSel.get(l2)]*ldMatOrg[snpSel.get(l1)][snpSel.get(l2)];
//			}
//		}
//		ldSum /= snpSel.size() * snpSel.size();
//		Logger.printUserLog("LD squres for all pair of markers is: " + ldSum);
//	}

	private void lineup()
	{
		//line up on the bim file
		boolean[] FileKeep = new boolean[archArgs.getMetaFile().length];
		Arrays.fill(FileKeep, true);

		gReader = new GWASReader(archArgs.getMetaFile(), FileKeep, archArgs.getKeys(), archArgs.isQT(), archArgs.isGZ(), archArgs.isChr(), archArgs.getChr());
		gReader.Start(false);

		HashMap<String, MetaStat> mStat = gReader.getMetaStat().get(0);

		effectSel = NewIt.newArrayList();
		snpSel = NewIt.newArrayList();

		int cntAmb = 0;
		int cntNotBiallele = 0;

		for(int i = 0; i < snpList.size(); i++)
		{
			SNP snp = snpList.get(i);
			if (SNPMatch.isAmbiguous(snp.getFirstAllele(), snp.getSecAllele())) 
			{
				cntAmb++;
				continue;
			}

			if (mStat.containsKey(snp.getName()))
			{
				MetaStat ms = mStat.get(snp.getName());

				if (!SNPMatch.IsBiallelic(snp.getFirstAllele(), snp.getSecAllele(), ms.getA1(), ms.getA2()))
				{
					cntNotBiallele++;
					continue;
				}

				if (SNPMatch.isAllelesMatchForTwoLoci(snp.getFirstAllele(), snp.getSecAllele(), ms.getA1(), ms.getA2()))
				{
					if (SNPMatch.isAlleleMatch(snp.getFirstAllele(), ms.getA1()))
					{
						effectSel.add((double) ms.getEffect());
					}
					else if (SNPMatch.isAlleleMatch(snp.getFirstAllele(), ms.getA2()))
					{
						effectSel.add((double) -1*ms.getEffect());
					}
				}
				else if (SNPMatch.isAllelesFlipMatchForTwoLoci(snp.getFirstAllele(), snp.getSecAllele(), ms.getA1(), ms.getA2()))
				{
					if (SNPMatch.isAlleleFlipMatch(snp.getFirstAllele(), ms.getA1()))
					{
						effectSel.add((double)ms.getEffect());
					}
					else if (SNPMatch.isAlleleFlipMatch(snp.getFirstAllele(), ms.getA2()))
					{
						effectSel.add((double)ms.getEffect());
					}
				}
				snpSel.add(i);
			}
		}
		Logger.printUserLog("Matched " + snpSel.size() + " markers.");

		if (cntAmb > 0)
		{
			Logger.printUserLog("Removed " + cntAmb + " parlindromic marker(s).");
		}

		if (cntNotBiallele > 0)
		{
			Logger.printUserLog("Removed " + cntNotBiallele + " biallelic marker(s).");
		}

		snpListAdj = NewIt.newArrayList();
		snpListAdjName = NewIt.newArrayList();
		for (int i = 0; i < snpSel.size(); i++)
		{
			snpListAdj.add(snpList.get(snpSel.get(i)));
			snpListAdjName.add(snpList.get(snpSel.get(i)).getName());
		}
	}

//	private void readLDscore()
//	{
//		BufferedReader reader = BufferedReader.openGZipFile(archArgs.getLDRFile(), "LD score (correlation) information");
//		ldMatOrg = new float[snpList.size()][snpList.size()];
//
//		String[] tokens = null;
//		long cnt=0;
//		while( (tokens = reader.readTokensAtLeast(4)) != null)
//		{
//			int i1 = Integer.parseInt(tokens[0])-1;
//			int i2 = Integer.parseInt(tokens[1])-1;
//			float score = Float.parseFloat(tokens[3]);
//			ldMatOrg[i1][i2] = ldMatOrg[i2][i1] = score;
//			cnt++;
//		}
//		reader.close();
//		Logger.printUserLog("Read " + cnt + " score from '" + archArgs.getLDRFile()+"'.");
//	}

	private void readSNPinfo()
	{
		BufferedReader reader = BufferedReader.openTextFile(archArgs.getBim(), "SNP information");
		snpList = NewIt.newArrayList();
		String[] tokens = null;
		while((tokens=reader.readTokensAtLeast(6))!= null)
		{
			String chr = tokens[0];
			String name = tokens[1];
			float dis = Float.parseFloat(tokens[2]);
			int pos = Integer.parseInt(tokens[3]);
			char refA = tokens[4].charAt(0);
			char altA = tokens[5].charAt(0);
			snpList.add(new SNP(chr, name, dis, pos, refA, altA));
		}
		reader.close();
		Logger.printUserLog("Read " + snpList.size() + " snps from '" + archArgs.getBim() + "'.");
	}

	private void readExtract()
	{
		BufferedReader reader = BufferedReader.openTextFile(archArgs.getExtractFile(), "Extract snps");
		snpExtract = NewIt.newArrayList();
		String[] tokens = null;
		while((tokens=reader.readTokensAtLeast(1))!= null)
		{
			String name = tokens[0];
			snpExtract.add(name);
		}
		reader.close();
		Logger.printUserLog("Read " + snpExtract.size() + " snps from '" + archArgs.getExtractFile() + "'.");
	}

	private ArchCommandArguments archArgs = null;
	private ArrayList<String> snpExtract = null;
	private ArrayList<SNP> snpList = null;
	
	private ArrayList<Integer> snpSel = null;
	private ArrayList<Double> effectSel = null;
	private GWASReader gReader;

	private ArrayList<SNP> snpListAdj = null;
	private ArrayList<String> snpListAdjName = null;
	
	private MapFile mapFile;
	private SampleFilter sf;
	private SumStatQC ssQC;
	private GenotypeMatrix gm;

}
