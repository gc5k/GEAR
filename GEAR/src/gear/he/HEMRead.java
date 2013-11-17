package gear.he;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.zip.GZIPInputStream;

import gear.CmdArgs;
import gear.ConstValues;
import gear.data.SubjectID;
import gear.he.covar.HeCov;
import gear.util.BinaryInputFile;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;

public class HEMRead
{
	protected boolean[] flag;
	HashMap<String, Integer> ID;
	protected String grmListFile;
	protected String grmFile;
	protected String grmID;
	protected String keepFile;
	protected String phenoFile;
	protected String output;
	protected int perm;
	protected boolean permFlag = false;

	protected boolean reverse;
	protected boolean k_button;
	protected double k;
	protected gear.HEType heType;

	protected double yyProd;

	protected double[][] XtX;
	protected double[] XtY;
	protected double[] y;

	protected ArrayList<BufferedReader> grmList;
	protected ArrayList<BinaryInputFile> binList;
	protected ArrayList<String> grmFileList;
	protected ArrayList<String> idFileList;
	
	protected boolean isSingleGrm;
	protected boolean isBinGrm = false;

	protected long nRec = 0;
	protected Lambda lambda;
	StringBuffer sb = new StringBuffer();

	public HEMRead()
	{
		if (CmdArgs.INSTANCE.getHEArgs().isSingleGrm())
		{
			grmFile = CmdArgs.INSTANCE.getHEArgs().getGrm();
			grmID = CmdArgs.INSTANCE.getHEArgs().getGrmId();
			existGrm();
			isSingleGrm = true;
		}
		else if (CmdArgs.INSTANCE.getHEArgs().isMultiGrm())
		{
			grmListFile = CmdArgs.INSTANCE.getHEArgs().getGrmList();
			existGrmList();
			isSingleGrm = false;
		}

		keepFile = CmdArgs.INSTANCE.keepFile;
		phenoFile = CmdArgs.INSTANCE.getHEArgs().getPheno();
		reverse = CmdArgs.INSTANCE.reverse;
		k_button = CmdArgs.INSTANCE.k_button;
		k = CmdArgs.INSTANCE.k;
		output = CmdArgs.INSTANCE.out;
		heType = CmdArgs.INSTANCE.getHEArgs().getType();
		permFlag = CmdArgs.INSTANCE.permFlag;
		perm = CmdArgs.INSTANCE.perm;

		HashMap<SubjectID, Integer> id2Idx = new HashMap<SubjectID, Integer>();
		
		readGrmIds(id2Idx);
		flag = new boolean[id2Idx.size()];
		Arrays.fill(flag, false);
		nRec = flag.length * (flag.length + 1) / 2;

		readPhenotypes(id2Idx);

		if (CmdArgs.INSTANCE.covar_file != null & CmdArgs.INSTANCE.qcovar_file != null) 
		{
			HeCov hecov = new HeCov(y, flag, id2Idx, CmdArgs.INSTANCE.qcovar_file, CmdArgs.INSTANCE.qcovar_num, CmdArgs.INSTANCE.covar_file, CmdArgs.INSTANCE.covar_num);
			hecov.generate_Res();
			y = hecov.getAdjustedPhe();
			flag = hecov.getFlag();
		}

		keepSpecifiedSubjects(id2Idx);

		int numAvailSubjects = 0;
		for (int i = 0; i < flag.length; i++)
		{
			if (flag[i])
			{
				++numAvailSubjects;
			}
		}
		
		Logger.printUserLog("Individuals matched phenotype and grm:" + numAvailSubjects);
		if (CmdArgs.INSTANCE.scale) 
		{
			Logger.printUserLog("standardizing the phenotype...");
			standardisePhenotypes(numAvailSubjects);
		}
	}


	private void existGrm()
	{
		grmList = NewIt.newArrayList();
		grmFileList = NewIt.newArrayList();
		idFileList = NewIt.newArrayList();
		binList = NewIt.newArrayList();

		if (CmdArgs.INSTANCE.getHEArgs().isGrmTxt())
		{
			FileUtil.exists(grmFile);
			FileUtil.exists(grmID);

			grmFileList.add(grmFile);
			idFileList.add(grmID);

			grmList.add(FileUtil.FileOpen(grmFile
					.toString()));
		}
		else if (CmdArgs.INSTANCE.getHEArgs().isGrm())
		{
			FileUtil.exists(grmFile);
			FileUtil.exists(grmID);

			grmFileList.add(grmFile);
			idFileList.add(grmID);

			FileInputStream fin = null;
			try
			{
				fin = new FileInputStream(grmFile);
			}
			catch (FileNotFoundException e1)
			{
				Logger.handleException(e1,
										"Error in opening '" + grmFile + "'.");
			}
			GZIPInputStream gzis = null;
			try
			{
				gzis = new GZIPInputStream(fin);
			}
			catch (IOException e1)
			{
				Logger.handleException(e1,
										"Error in opening gz file '" + grmFile + "'.");
			}
			InputStreamReader xover = new InputStreamReader(gzis);

			BufferedReader grmFile = new BufferedReader(xover);
			grmList.add(grmFile);
		}
		else
		{
			FileUtil.exists(grmFile);
			FileUtil.exists(grmID);

			grmFileList.add(grmFile);
			idFileList.add(grmID);

			BinaryInputFile grmBin = new BinaryInputFile(grmFile, "GRM");
			grmBin.setLittleEndian(true);
			binList.add(grmBin);
		}
		if (CmdArgs.INSTANCE.getHEArgs().isGrm() || CmdArgs.INSTANCE
				.getHEArgs().isGrmTxt())
		{
			if (grmList.size() == 0)
			{
				Logger.printUserError("Empty '" + CmdArgs.INSTANCE.getHEArgs()
						.getGrmList() + "'.");
				System.exit(1);
			}
		}
		if (CmdArgs.INSTANCE.getHEArgs().isGrmBinary())
		{
			if (binList.size() == 0)
			{
				Logger.printUserError("Empty '" + CmdArgs.INSTANCE.getHEArgs()
						.getGrmList() + "'.");
				System.exit(1);
			}
			isBinGrm = true;
		}
	}
	
	public boolean isCaseControl()
	{
		return isCaseCtrl;
	}
	
	public double getCaseValue()
	{
		return caseVal;
	}
	
	public double getControlValue()
	{
		return ctrlVal;
	}
	
	public int getNumberOfCases()
	{
		return numCases;
	}
	
	public int getNumberOfControls()
	{
		return numCtrls;
	}

	private void existGrmList()
	{
		grmList = NewIt.newArrayList();
		grmFileList = NewIt.newArrayList();
		idFileList = NewIt.newArrayList();
		binList = NewIt.newArrayList();

		BufferedReader listReader = FileUtil.FileOpen(grmListFile);
		String line = null;
		try
		{
			while ((line = listReader.readLine()) != null)
			{
				String root = line.trim();
				StringBuilder tsb = new StringBuilder(root);

				if (CmdArgs.INSTANCE.getHEArgs().isGrmTxtList())
				{
					FileUtil.exists(tsb.append(".grm.txt").toString());
					tsb = new StringBuilder(root);
					FileUtil.exists(tsb.append(".grm.id").toString());

					tsb = new StringBuilder(root);
					grmFileList.add(tsb.append(".grm.txt").toString());
					tsb = new StringBuilder(root);
					idFileList.add(tsb.append(".grm.id").toString());

					tsb = new StringBuilder(root);
					grmList.add(FileUtil.FileOpen(tsb.append(".grm.txt")
							.toString()));
				}
				else if (CmdArgs.INSTANCE.getHEArgs().isGrmList())
				{
					FileUtil.exists(tsb.append(".grm.gz").toString());
					tsb = new StringBuilder(root);
					FileUtil.exists(tsb.append(".grm.id").toString());

					tsb = new StringBuilder(root);
					grmFileList.add(tsb.append(".grm.gz").toString());
					tsb = new StringBuilder(root);
					idFileList.add(tsb.append(".grm.id").toString());

					FileInputStream fin = null;
					try
					{
						tsb = new StringBuilder(root);
						fin = new FileInputStream(tsb.append(".grm.gz")
								.toString());
					}
					catch (FileNotFoundException e1)
					{
						Logger.handleException(	e1,
												"Error in opening '" + line
														.trim() + "' in " + CmdArgs.INSTANCE
														.getHEArgs()
														.getGrmList() + "'.");
					}
					GZIPInputStream gzis = null;
					try
					{
						gzis = new GZIPInputStream(fin);
					}
					catch (IOException e1)
					{
						Logger.handleException(	e1,
												"Error in opening gz file '" + line
														.trim() + "' in " + CmdArgs.INSTANCE
														.getHEArgs()
														.getGrmList() + "'.");
					}
					InputStreamReader xover = new InputStreamReader(gzis);

					BufferedReader grmFile = new BufferedReader(xover);
					grmList.add(grmFile);
				}
				else
				{
					tsb = new StringBuilder(root);
					FileUtil.exists(tsb.append(".grm.bin").toString());
					tsb = new StringBuilder(root);
					FileUtil.exists(tsb.append(".grm.id").toString());

					tsb = new StringBuilder(root);
					grmFileList.add(tsb.append(".grm.bin").toString());
					tsb = new StringBuilder(root);
					idFileList.add(tsb.append(".grm.id").toString());

					BinaryInputFile grmBin = new BinaryInputFile(root + ".grm.bin", "GRM");
					grmBin.setLittleEndian(true);
					binList.add(grmBin);
				}
			}
		}
		catch (IOException e)
		{
			Logger.handleException(	e,
									"An exception occurred when reading the GRM file '" + CmdArgs.INSTANCE
											.getHEArgs().getGrmList() + "'.");
		}
		if (CmdArgs.INSTANCE.getHEArgs().isGrmList() || CmdArgs.INSTANCE
				.getHEArgs().isGrmTxtList())
		{
			if (grmList.size() == 0)
			{
				Logger.printUserError("Empty '" + CmdArgs.INSTANCE.getHEArgs()
						.getGrmList() + "'.");
				System.exit(1);
			}
		}
		if (CmdArgs.INSTANCE.getHEArgs().isGrmBinaryList())
		{
			if (binList.size() == 0)
			{
				Logger.printUserError("Empty '" + CmdArgs.INSTANCE.getHEArgs()
						.getGrmList() + "'.");
				System.exit(1);
			}
			isBinGrm = true;
		}
		grmID = idFileList.get(0);
	}
	
	private void readGrmIds(HashMap<SubjectID, Integer> id2Idx)
	{

		gear.util.BufferedReader reader = gear.util.BufferedReader.openTextFile(grmID, "GRM-ID");
		int idx = 0;
		String[] tokens;
		while ((tokens = reader.readTokens(2)) != null)
		{
			id2Idx.put(new SubjectID(tokens[0], tokens[1]), idx++);
		}
		Logger.printUserLog("individuals in grm id file: " + id2Idx.size());
		reader.close();

	}
	
	private void readPhenotypes(HashMap<SubjectID, Integer> id2Idx)
	{
		gear.util.BufferedReader reader = gear.util.BufferedReader.openTextFile(phenoFile, "phenotype");
		
		y = new double[id2Idx.size()];
		
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
		
		checkCaseControl();
	}
	
	private void checkCaseControl()
	{
		double val1 = 0.0, val2 = 0.0;
		int numVal1 = 0, numVal2 = 0;
		
		for (int ii = 0; ii < y.length; ++ii)
		{
			if (flag[ii])
			{
				if (numVal1 == 0)
				{
					val1 = y[ii];
					++numVal1;
				}
				else if (Math.abs(val1 - y[ii]) <= ConstValues.EPSILON)
				{
					++numVal1;
				}
				else if (numVal2 == 0)
				{
					val2 = y[ii];
					++numVal2;
				}
				else if (Math.abs(val2 - y[ii]) <= ConstValues.EPSILON)
				{
					++numVal2;
				}
				else
				{
					isCaseCtrl = false;
					return;
				}
			}
		}
		
		isCaseCtrl = true;
		
		if (val1 > val2)
		{
			caseVal = val1;
			ctrlVal = val2;
			numCases = numVal1;
			numCtrls = numVal2;
		}
		else
		{
			caseVal = val2;
			ctrlVal = val1;
			numCases = numVal2;
			numCtrls = numVal1;
		}
	}
	
	private void standardisePhenotypes(int numAvailSubjects)
	{
		Logger.printUserLog("Standardising phentoypes.");

		double ss = 0.0, ssx = 0.0;
		for (int i = 0; i < flag.length; i++)
		{
			if (flag[i])
			{
				ss += y[i];
				ssx += y[i] * y[i];
			}
		}
		ss /= numAvailSubjects;
		
		double sd = Math.sqrt((ssx - numAvailSubjects * ss * ss) / (numAvailSubjects - 1));


		for (int i = 0; i < flag.length; i++)
		{
			if (flag[i])
			{
				y[i] = (y[i] - ss) / sd;
			}
		}
	}
	
	private void keepSpecifiedSubjects(HashMap<SubjectID, Integer> id2Idx)
	{
		if (keepFile != null)
		{
			gear.util.BufferedReader reader = gear.util.BufferedReader.openTextFile(keepFile, "keep-individual");
			boolean[] ff = new boolean[flag.length];
			Arrays.fill(ff, false);

			String[] tokens = null;
			while ((tokens = reader.readTokens(2)) != null)
			{
				Integer idx = id2Idx.get(new SubjectID(/*famID*/tokens[0], /*indID*/tokens[1]));
				if (idx != null)
				{
					ff[idx] = true;
				}
			}
			reader.close();

			for (int i = 0; i < ff.length; i++)
			{
				flag[i] &= ff[i];
			}
		}
	}
	
	// Case-Control Cached Data
	private boolean isCaseCtrl;
	private double caseVal, ctrlVal;
	private int numCases, numCtrls;
}
