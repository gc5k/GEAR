package gear.bluppca;

import gear.CommandArguments;
import gear.CommandImpl;
import gear.ConstValues;
import gear.family.pedigree.file.MapFile;
import gear.family.pedigree.file.SNP;
import gear.family.plink.PLINKParser;
import gear.family.popstat.GenotypeMatrix;
import gear.family.qc.rowqc.SampleFilter;
import gear.he.SubjectID;
import gear.sumstat.qc.rowqc.SumStatQC;
import gear.util.BinaryInputFile;
import gear.util.BufferedReader;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.pop.PopStat;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;

public class BlupPCACommandImpl extends CommandImpl
{
	@Override
	public void execute(CommandArguments cmdArgs)
	{
		BlupPCACommandArguments blupArgs = (BlupPCACommandArguments)cmdArgs;

		//read grm
		id2Idx = new HashMap<SubjectID, Integer>();
		readGRM_IDs(blupArgs.getGRM_ID(), id2Idx);
		readGRM(blupArgs);

		flag = new boolean[id2Idx.size()];
		Arrays.fill(flag, false);

		readPhenotypes(blupArgs.getPhenotypeFile(), id2Idx);

		PLINKParser pp = PLINKParser.parse(blupArgs);
		sf = new SampleFilter(pp.getPedigreeData(), pp.getMapData());
		ssQC = new SumStatQC(pp.getPedigreeData(), pp.getMapData(), sf);
		mapFile = ssQC.getMapFile();
		gm = new GenotypeMatrix(ssQC.getSample());

		PopStat.Imputation(gm);

		double[][] genoMat = new double[gm.getNumIndivdial()][gm.getNumMarker()];
		for(int i = 0; i < genoMat.length; i++)
		{
			for(int j = 0; j < genoMat[i].length; j++)
			{
				genoMat[i][j] = gm.getAdditiveScoreOnFirstAllele(i, j);
			}
		}

		double[][] blupPC = new double[gm.getNumMarker()][phe[0].length];

		RealMatrix grm = new Array2DRowRealMatrix(A);
		RealMatrix grm_Inv = (new LUDecompositionImpl(grm)).getSolver().getInverse();
		
		Logger.printUserLog("Revving up the BLUP machine...");
		RealMatrix tmp = (new Array2DRowRealMatrix(genoMat)).transpose().multiply(grm_Inv);

		for(int i = 0; i < phe[0].length; i++)
		{
			Logger.printUserLog("Calculating blup vector[" + (i+1) + "].");

			double[] Y = new double[phe.length];
			for(int j = 0; j < phe.length; j++)
			{
				Y[j] = phe[j][i];
			}
			RealMatrix B = tmp.multiply(new Array2DRowRealMatrix(Y));

			Logger.printUserLog("Rescaling the snp effects...");
			for(int j = 0; j < B.getRowDimension(); j++)
			{
				blupPC[j][i] = B.getEntry(j, 0);
			}
		}

		PrintStream predictorFile = FileUtil.CreatePrintStream(blupArgs.getOutRoot() + ".blup");

		// Title Line
		ArrayList<SNP> snpList = mapFile.getMarkerList();

		predictorFile.print("SNP\tRefAllele");
		for(int i = 0; i < phe[0].length; i++)
		{
			if (i == (phe[0].length - 1)) 
			{
				predictorFile.println("\tBLUP" + (i+1));				
			}
			else
			{
				predictorFile.print("\tBLUP" + (i+1));
			}
		}

		for(int i = 0; i < gm.getNumMarker(); i++)
		{
			SNP snp = snpList.get(i);
			predictorFile.print(snp.getName() + "\t" + snp.getFirstAllele() + "\t");
			for(int j = 0; j < blupPC[i].length; j++)
			{
				if (j == (blupPC[i].length - 1))
				{
					predictorFile.println(blupPC[i][j]);
				}
				else
				{
					predictorFile.print(blupPC[i][j]+"\t");
				}
			}
		}
		predictorFile.close();
	}

	private void readGRM_IDs(String fileName, HashMap<SubjectID, Integer> id2Idx)
	{
		BufferedReader reader = BufferedReader.openTextFile(fileName, "GRM-ID");
		int idx = 0;
		String[] tokens;
		while ((tokens = reader.readTokens(2)) != null)
		{
			id2Idx.put(new SubjectID(/*famID*/tokens[0], /*indID*/tokens[1]), idx++);
		}
		Logger.printUserLog("individuals in grm id file: " + id2Idx.size());
		reader.close();
	}
	
	private void readGRM(BlupPCACommandArguments blupArgs)
	{
		if (blupArgs.getGRMBin() != null)
		{
			readGRMBin(blupArgs.getGRMBin());
		}
		else
		{
			BufferedReader reader = blupArgs.getGRMText() == null ?
					BufferedReader.openGZipFile(blupArgs.getGRM_GZ(), "GRM (GZip)") :
					BufferedReader.openTextFile(blupArgs.getGRMText(), "GRM");
			readGRM(reader);
		}
	}	

	private void readGRMBin(String fileName)
	{
		BinaryInputFile grmBin = new BinaryInputFile(fileName, "GRM (Binary)", /*littleEndian*/true);
		A = new double[id2Idx.size()][id2Idx.size()];
		Logger.printUserLog("Constructing A matrix: a " + id2Idx.size() + " X " + id2Idx.size() + " matrix.");
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

	private void readGRM(BufferedReader reader)
	{
		A = new double[id2Idx.size()][id2Idx.size()];
		String[] tokens = null;
		for (int i = 0; i < A.length; i++)
		{
			for (int j = 0; j <= i; j++)
			{
				if ((tokens = reader.readTokens(4)) != null) 
				{
					A[i][j] = A[j][i] = Double.parseDouble(tokens[3]);
				}
			}
		}
		reader.close();
	}

	private void readPhenotypes(String fileName, HashMap<SubjectID, Integer> id2Idx)
	{
		gear.util.BufferedReader reader = gear.util.BufferedReader.openTextFile(fileName, "phenotype");
				
		@SuppressWarnings("unchecked")
		HashMap<SubjectID, Integer> subjectsUnread = (HashMap<SubjectID, Integer>)id2Idx.clone();

		HashSet<SubjectID> subjectsRead = new HashSet<SubjectID>();

		String[] tokens = null;

		int c = 0;
		while ((tokens = reader.readTokens()) != null)
		{

			if (tokens.length < 3)
			{
				reader.errorPreviousLine("There should be at least " + 3 + " columns.");
			}
			if (c == 0) 
			{
				phe = new double[id2Idx.size()][tokens.length - 2];
				c++;
				for( int i = 0; i < id2Idx.size(); i++)
				{
					Arrays.fill(phe[i], -9);
				}
			}

			SubjectID subID = new SubjectID(/*famID*/tokens[0], /*indID*/tokens[1]);

			int ii = 0;
			if (subjectsUnread.containsKey(subID))
			{
				ii = subjectsUnread.get(subID);
				boolean f = true;
				String pheValStr = null;
				try
				{
					for( int i = 0; i < phe[ii].length; i++)
					{
						pheValStr = tokens[2 + i];
						if (ConstValues.isNA(pheValStr))
						{
							f = false;
							break;
						}
						phe[ii][i] = Double.parseDouble(pheValStr);							
					}
				}
				catch (NumberFormatException e)
				{
					reader.errorPreviousLine("'" + pheValStr + "' is not a valid phenotype value. It should be a floating point number.");
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

	private boolean[] flag;
	private double[][] phe;
	private double[][] A;

	private HashMap<SubjectID, Integer> id2Idx;

	private MapFile mapFile;
	private SampleFilter sf;
	private SumStatQC ssQC;
	private GenotypeMatrix gm;
}
