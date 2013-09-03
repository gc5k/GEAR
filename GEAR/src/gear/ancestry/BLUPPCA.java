package gear.ancestry;

import gear.CmdArgs;
import gear.ConstValues;
import gear.family.pedigree.file.MapFile;
import gear.family.pedigree.file.SNP;
import gear.family.plink.PLINKBinaryParser;
import gear.family.plink.PLINKParser;
import gear.family.popstat.GenotypeMatrix;
import gear.family.qc.rowqc.SampleFilter;
import gear.he.SubjectID;
import gear.sumstat.qc.rowqc.SumStatQC;
import gear.util.BinaryInputFile;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.pop.PopStat;
import gear.util.stat.PCA.NCGStatUtils;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;

public class BLUPPCA
{
	
	private boolean[] flag;
	private String phenoFile;
	private double[] y;
	private double[][] A;

	private HashMap<SubjectID, Integer> id2Idx;

	private MapFile mapFile;
	private SampleFilter sf;
	private SumStatQC ssQC;
	private GenotypeMatrix gm;

	private double[][] freq;
	
	public BLUPPCA()
	{

		//read grm
		id2Idx = new HashMap<SubjectID, Integer>();
		readGrmIds(id2Idx);
		readGRM();

		flag = new boolean[id2Idx.size()];
		Arrays.fill(flag, false);

		//read pheno
		phenoFile = CmdArgs.INSTANCE.getHEArgs().getPheno();
		readPhenotypes(id2Idx);

		PLINKParser pp = null;
		if (CmdArgs.INSTANCE.getFileArgs().isSet())
		{
			pp = new PLINKParser(CmdArgs.INSTANCE.getFileArgs()
					.getPed(), CmdArgs.INSTANCE.getFileArgs()
					.getMap());
		}
		else if (CmdArgs.INSTANCE.getBFileArgs(0).isSet())
		{
			pp = new PLINKBinaryParser(CmdArgs.INSTANCE.getBFileArgs(0)
					.getBed(), CmdArgs.INSTANCE.getBFileArgs(0)
					.getBim(), CmdArgs.INSTANCE.getBFileArgs(0)
					.getFam());
		}
		else
		{
			Logger.printUserError("No input files.");
			System.exit(1);
		}
		pp.Parse();

		sf = new SampleFilter(pp.getPedigreeData(),
				pp.getMapData());
		ssQC = new SumStatQC(pp.getPedigreeData(), pp.getMapData(),
				sf);
		mapFile = ssQC.getMapFile();
		gm = new GenotypeMatrix(ssQC.getSample());
		//impute
		PopStat.Imputation(gm);
		//calculate freq
		freq = PopStat.calAlleleFrequency(gm, gm.getNumMarker());
	}

	public void BLUPit()
	{
		RealMatrix grm = new Array2DRowRealMatrix(A);
		RealMatrix grm_Inv = (new LUDecompositionImpl(grm)).getSolver().getInverse();

		double[][] genoMat = new double[gm.getNumIndivdial()][gm.getNumMarker()];
		for(int i = 0; i < genoMat.length; i++)
		{
			for(int j = 0; j < genoMat[i].length; j++)
			{
				genoMat[i][j] = gm.getAdditiveScore(i, j);
			}
		}

		Logger.printUserLog("Standardizing genotypes...");
		NCGStatUtils.standardize(genoMat, false);

		Logger.printUserLog("Calculating blup pca...");
		RealMatrix tmp = (new Array2DRowRealMatrix(genoMat)).transpose().multiply(grm_Inv);
		RealMatrix B = tmp.multiply(new Array2DRowRealMatrix(y));

		Logger.printUserLog("Rescaling the snp effects...");
		for(int i = 0; i < B.getRowDimension(); i++)
		{
			double b = B.getEntry(i, 0);
			double sd = Math.sqrt(freq[i][0] * (1-freq[i][0]) * 2);
			double s = b / sd / B.getRowDimension();
			B.setEntry(i, 0, s);
		}

		StringBuilder fsb = new StringBuilder();
		fsb.append(CmdArgs.INSTANCE.out);
		fsb.append(".blup");
		PrintStream predictorFile = FileUtil.CreatePrintStream(fsb.toString());

		// Title Line
		ArrayList<SNP> snpList = mapFile.getMarkerList();

		predictorFile.println("SNP\tA1\tBLUP");
		for(int i = 0; i < B.getRowDimension(); i++)
		{
			SNP snp = snpList.get(i);
			predictorFile.println(snp.getName() + "\t" + snp.getFirstAllele() + "\t" + B.getEntry(i, 0));
		}
		predictorFile.close();

	}

	private void readGrmIds(HashMap<SubjectID, Integer> id2Idx)
	{

		gear.util.BufferedReader reader = gear.util.BufferedReader.openTextFile(CmdArgs.INSTANCE.getHEArgs().getGrmId(), "GRM-ID");
		int idx = 0;
		String[] tokens;
		while ((tokens = reader.readTokens(2)) != null)
		{
			id2Idx.put(new SubjectID(tokens[0], tokens[1]), idx++);
		}
		Logger.printUserLog("individuals in grm id file: " + id2Idx.size());
		reader.close();
	}

	private void readGRM()
	{
		BinaryInputFile grmBin = new BinaryInputFile(CmdArgs.INSTANCE.getHEArgs().getGrm(), "GRM");
		grmBin.setLittleEndian(true);
		A = new double[id2Idx.size()][id2Idx.size()];
		Logger.printUserLog("making A matrix");
		for (int i = 0; i < A.length; i++) 
		{
			for (int j = 0; j <= i; j++)
			{
				if (grmBin.available() >= ConstValues.FLOAT_SIZE)
				{
					A[i][j] = A[j][i] = grmBin.readFloat();
				}
			}
		}
	}

	private void readPhenotypes(HashMap<SubjectID, Integer> id2Idx)
	{
		Logger.printUserLog("reading phentoypes from '" + phenoFile + "'");
		gear.util.BufferedReader reader = gear.util.BufferedReader.openTextFile(phenoFile, "phenotype");
		
		y = new double[id2Idx.size()];
		Arrays.fill(y, -9);
		
		@SuppressWarnings("unchecked")
		HashMap<SubjectID, Integer> subjectsUnread = (HashMap<SubjectID, Integer>)id2Idx.clone();

		HashSet<SubjectID> subjectsRead = new HashSet<SubjectID>();

		int tarTraitIdx = CmdArgs.INSTANCE.getHEArgs().getTargetTraitOptionValue() - 1;
		int minNumCols = 2 + tarTraitIdx + 1;

		String[] tokens = null;

		while ((tokens = reader.readTokens()) != null)
		{
			if (tokens.length < minNumCols)
			{
				reader.errorPreviousLine("There should be at least " + minNumCols + " columns.");
			}
			
			SubjectID subID = new SubjectID(/*famID*/tokens[0], /*indID*/tokens[1]);

			int ii = 0;
			if (subjectsUnread.containsKey(subID))
			{
				ii = subjectsUnread.get(subID);
				boolean f = true;
				String pheValStr = tokens[2 + tarTraitIdx];
				if (ConstValues.isNA(pheValStr))
				{
					f = false;
					break;
				}
				else
				{
					try
					{
						y[ii] = Double.parseDouble(pheValStr);
					}
					catch (NumberFormatException e)
					{
						reader.errorPreviousLine("'" + pheValStr + "' is not a valid phenotype value. It should be a floating point number.");
					}
				}
				flag[ii] = f;
				
				subjectsUnread.remove(subID);
				subjectsRead.add(subID);
			}
			else if (subjectsRead.contains(subID))
			{
				reader.errorPreviousLine("Individual " + subID + " is repeated.");
			}
			else
			{
				flag[ii] = false;
//				reader.errorPreviousLine("Individual " + subID + " appears in the phenotype file but not in the grm id file(s).");
			}
		}
		reader.close();

		if (!subjectsUnread.isEmpty())
		{
			String msg = "";
			msg += subjectsUnread.size() + " individual(s) (e.g. " + subjectsUnread.keySet().iterator().next();
			msg += ") appear in the grm id file(s) but not in the phenotype file";
			Logger.printUserError(msg);
		}
	}
}
